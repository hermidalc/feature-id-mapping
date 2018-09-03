#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use File::Basename qw(fileparse);
use Getopt::Long qw(:config auto_help auto_version);
use Pod::Usage qw(pod2usage);
use Sort::Key qw(nsort);
use Sort::Key::Natural qw(natsort);
use Storable qw(lock_retrieve);
use Data::Dumper;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Terse = 1;
$Data::Dumper::Deepcopy = 1;

my $src_id_type = '';
my $map_id_type = 'hgnc_symbol';
my $map_id_file = '';
my $debug = 0;
my $verbose = 0;
GetOptions(
    'src-id-type:s'     => \$src_id_type,
    'map-id-type:s'     => \$map_id_type,
    'map-id-file:s'     => \$map_id_file,
    'debug'             => \$debug,
    'verbose'           => \$verbose,
) || pod2usage(-verbose => 0);
pod2usage(-message => 'Missing source ID file param') unless @ARGV and scalar(@ARGV) == 1;
pod2usage(-message => "'$ARGV[0]' is not a valid source ID file path") unless -f $ARGV[0];
pod2usage(-message => 'Missing required --map-id-file=<file>') if $map_id_type ne 'hgnc_symbol' and !$map_id_file;
# load reference data
my $organism_name = 'Homo sapiens';
my $CTK_DATA_ID_MAPPING_GENE_SYMBOL_SUFFIX = '_gene_symbols';
my $gene_info_hashref = lock_retrieve("$FindBin::Bin/gene_info.pls");
my $gene_history_hashref = lock_retrieve("$FindBin::Bin/gene_history.pls");
my $gene_history_by_symbol_hashref = lock_retrieve("$FindBin::Bin/gene_history_by_symbol.pls");
my $ensembl2gene_ids_map = lock_retrieve("$FindBin::Bin/ensembl2gene_ids.pls");
my $gene2ensembl_ids_map = lock_retrieve("$FindBin::Bin/gene2ensembl_ids.pls");
(my $organism_file_basename = $organism_name) =~ s/\s+/_/g;
my $symbol2gene_ids_map =
    lock_retrieve("$FindBin::Bin/${organism_file_basename}${CTK_DATA_ID_MAPPING_GENE_SYMBOL_SUFFIX}.map.pls");
my $gene_id2symbols_map =
    lock_retrieve("$FindBin::Bin/${organism_file_basename}${CTK_DATA_ID_MAPPING_GENE_SYMBOL_SUFFIX}.revmap.pls");
my $gene_symbol2id_bestmap =
    lock_retrieve("$FindBin::Bin/${organism_file_basename}${CTK_DATA_ID_MAPPING_GENE_SYMBOL_SUFFIX}.bestmap.pls");
my $find_entrez_map = $map_id_type eq 'hgnc_symbol' ? 1 : 0;
my $out_file_suffix = 'hgnc_symbols';
my (@ids, %id_map, %mapped_ids, %unmapped_ids, %custom_map_ids, %ambig_map);
# load custom ids
if ($map_id_type ne 'hgnc_symbol') {
    ($out_file_suffix, undef, undef) = fileparse($map_id_file, qr/\.[^.]*/);
    open(my $map_fh, '<', $map_id_file);
    while (my $id = <$map_fh>) {
        if ($map_id_type eq 'custom_entrez') {
            $id =~ s/\D//g;
        }
        elsif ($map_id_type eq 'custom_ensembl') {
            $id =~ s/\s+//g;
        }
        elsif ($map_id_type eq 'custom_symbol') {
            $id =~ s/\s+//g;
            if (
                !exists($gene_symbol2id_bestmap->{$id}) or
                exists($gene_symbol2id_bestmap->{$id}->{ambig_gene_map})
            ) {
                # print "$id\n";
            }
        }
        $custom_map_ids{$id}++;
    }
    close($map_fh);
}
my ($file_basename, undef, $file_ext) = fileparse($ARGV[0], qr/\.[^.]*/);
open(my $in_fh, '<', $ARGV[0]) or die "$!";
while (my $line = <$in_fh>) {
    my ($id) = split /\t/, $line;
    if ($src_id_type eq 'entrez') {
        $id =~ s/\D//g;
    }
    elsif ($src_id_type eq 'ensembl') {
        $id =~ s/\s+//g;
    }
    elsif ($src_id_type eq 'symbol') {
        $id =~ s/\s+//g;
    }
    push @ids, $id;
}
close($in_fh);
# entrez
if ($src_id_type eq 'entrez') {
    for my $entrez_gene_id (@ids) {
        my $mapped_symbol = find_entrez_gene_id_symbol($entrez_gene_id);
        if (defined($mapped_symbol)) {
            $id_map{$entrez_gene_id} = $mapped_symbol;
            $mapped_ids{$mapped_symbol}++;
        }
        else {
            $unmapped_ids{$entrez_gene_id}++;
        }
    }
}
# ensembl
elsif ($src_id_type eq 'ensembl') {
    my $ensembl2gene_93_map;
    open(my $ensembl_fh, '<', "$FindBin::Bin/ensembl2gene_93_map.txt");
    while (my $line = <$ensembl_fh>) {
        chomp($line);
        my ($ensembl_stable_gene_id, $entrez_gene_id, $symbol) = split(/\t/, $line, 3);
        if (
            (defined($entrez_gene_id) and $entrez_gene_id ne '') or
            (defined($symbol) and $symbol ne '')
        ) {
            $entrez_gene_id = '' if not defined $entrez_gene_id;
            $symbol = '' if not defined $symbol;
            $ensembl2gene_93_map->{$ensembl_stable_gene_id}->{$symbol}->{$entrez_gene_id}++;
        }
    }
    close($ensembl_fh);
    # for my $ensembl_stable_gene_id (natsort keys %{$ensembl2gene_93_map}) {
    #     if (scalar(keys %{$ensembl2gene_93_map->{$ensembl_stable_gene_id}}) > 1) {
    #         print "'$ensembl_stable_gene_id' => ", Dumper($ensembl2gene_93_map->{$ensembl_stable_gene_id}), "\n";
    #     }
    # }
    my $ensembl2gene_grch37_map;
    open($ensembl_fh, '<', "$FindBin::Bin/ensembl2gene_grch37_map.txt");
    while (my $line = <$ensembl_fh>) {
        chomp($line);
        my ($ensembl_stable_gene_id, $entrez_gene_id, $symbol) = split(/\t/, $line, 3);
        if (
            (defined($entrez_gene_id) and $entrez_gene_id ne '') or
            (defined($symbol) and $symbol ne '')
        ) {
            $entrez_gene_id = '' if not defined $entrez_gene_id;
            $symbol = '' if not defined $symbol;
            $ensembl2gene_grch37_map->{$ensembl_stable_gene_id}->{$symbol}->{$entrez_gene_id}++;
        }
    }
    close($ensembl_fh);
    # for my $ensembl_stable_gene_id (natsort keys %{$ensembl2gene_grch37_map}) {
    #     if (scalar(keys %{$ensembl2gene_grch37_map->{$ensembl_stable_gene_id}}) > 1) {
    #         print "'$ensembl_stable_gene_id' => ", Dumper($ensembl2gene_grch37_map->{$ensembl_stable_gene_id}), "\n";
    #     }
    # }
    for my $ensembl_gene_id (@ids) {
        (my $ensembl_stable_gene_id = $ensembl_gene_id) =~ s/\.\d+$//;
        if (exists($ensembl2gene_93_map->{$ensembl_stable_gene_id})) {
            if (scalar(keys %{$ensembl2gene_93_map->{$ensembl_stable_gene_id}}) == 1) {
                my ($mapped_symbol) = keys %{$ensembl2gene_93_map->{$ensembl_stable_gene_id}};
                if ($mapped_symbol ne '') {
                    if (
                        $map_id_type eq 'hgnc_symbol' or
                        ($map_id_type eq 'custom_symbol' and exists($custom_map_ids{$mapped_symbol}))
                    ) {
                        $id_map{$ensembl_gene_id} = $mapped_symbol;
                        $mapped_ids{$mapped_symbol}++;
                    }
                }
                elsif (scalar(keys %{$ensembl2gene_93_map->{$ensembl_stable_gene_id}->{$mapped_symbol}}) == 1) {
                    my ($entrez_gene_id) = keys %{$ensembl2gene_93_map->{$ensembl_stable_gene_id}->{$mapped_symbol}};
                    if ($entrez_gene_id ne '') {
                        my $mapped_symbol = find_entrez_gene_id_symbol($entrez_gene_id);
                        if (defined($mapped_symbol) and (
                            $map_id_type eq 'hgnc_symbol' or
                            ($map_id_type eq 'custom_symbol' and exists($custom_map_ids{$mapped_symbol}))
                        )) {
                            $id_map{$ensembl_gene_id} = $mapped_symbol;
                            $mapped_ids{$mapped_symbol}++;
                        }
                    }
                }
                else {
                    push @{$ambig_map{$ensembl_gene_id}},
                        keys %{$ensembl2gene_93_map->{$ensembl_stable_gene_id}->{$mapped_symbol}};
                }
            }
            else {
                push @{$ambig_map{$ensembl_gene_id}},
                    keys %{$ensembl2gene_93_map->{$ensembl_stable_gene_id}};
            }
        }
        if (!exists($id_map{$ensembl_gene_id}) and exists($ensembl2gene_ids_map->{$ensembl_stable_gene_id})) {
            if (scalar(keys %{$ensembl2gene_ids_map->{$ensembl_stable_gene_id}}) == 1) {
                my ($entrez_gene_id) = keys %{$ensembl2gene_ids_map->{$ensembl_stable_gene_id}};
                my $mapped_symbol = find_entrez_gene_id_symbol($entrez_gene_id);
                if (defined($mapped_symbol) and (
                    $map_id_type eq 'hgnc_symbol' or
                    ($map_id_type eq 'custom_symbol' and exists($custom_map_ids{$mapped_symbol}))
                )) {
                    $id_map{$ensembl_gene_id} = $mapped_symbol;
                    $mapped_ids{$mapped_symbol}++;
                }
                delete $ambig_map{$ensembl_gene_id};
            }
            else {
                for my $entrez_gene_id (nsort keys %{$ensembl2gene_ids_map->{$ensembl_stable_gene_id}}) {
                    if (scalar(keys %{$gene2ensembl_ids_map->{$entrez_gene_id}}) == 1) {
                        my $mapped_symbol = find_entrez_gene_id_symbol($entrez_gene_id);
                        if (defined($mapped_symbol) and (
                            $map_id_type eq 'hgnc_symbol' or
                            ($map_id_type eq 'custom_symbol' and exists($custom_map_ids{$mapped_symbol}))
                        )) {
                            $id_map{$ensembl_gene_id} = $mapped_symbol;
                            $mapped_ids{$mapped_symbol}++;
                            delete $ambig_map{$ensembl_gene_id};
                            last;
                        }
                    }
                }
            }
        }
        if (!exists($id_map{$ensembl_gene_id}) and exists($ensembl2gene_grch37_map->{$ensembl_stable_gene_id})) {
            if (scalar(keys %{$ensembl2gene_grch37_map->{$ensembl_stable_gene_id}}) == 1) {
                my ($mapped_symbol) = keys %{$ensembl2gene_grch37_map->{$ensembl_stable_gene_id}};
                if ($mapped_symbol ne '') {
                    if ($map_id_type eq 'hgnc_symbol' or (
                        $map_id_type eq 'custom_symbol' and
                        exists($custom_map_ids{$mapped_symbol}) and
                        !exists($mapped_ids{$mapped_symbol})
                    )) {
                        $id_map{$ensembl_gene_id} = $mapped_symbol;
                        $mapped_ids{$mapped_symbol}++;
                    }
                }
                elsif (scalar(keys %{$ensembl2gene_grch37_map->{$ensembl_stable_gene_id}->{$mapped_symbol}}) == 1) {
                    my ($entrez_gene_id) = keys %{$ensembl2gene_grch37_map->{$ensembl_stable_gene_id}->{$mapped_symbol}};
                    if ($entrez_gene_id ne '') {
                        my $mapped_symbol = find_entrez_gene_id_symbol($entrez_gene_id);
                        if (defined($mapped_symbol) and ($map_id_type eq 'hgnc_symbol' or (
                            $map_id_type eq 'custom_symbol' and
                            exists($custom_map_ids{$mapped_symbol}) and
                            !exists($mapped_ids{$mapped_symbol})
                        ))) {
                            $id_map{$ensembl_gene_id} = $mapped_symbol;
                            $mapped_ids{$mapped_symbol}++;
                        }
                    }
                }
                elsif (!exists($ambig_map{$ensembl_gene_id})) {
                    push @{$ambig_map{$ensembl_gene_id}},
                        keys %{$ensembl2gene_grch37_map->{$ensembl_stable_gene_id}->{$mapped_symbol}};
                }
            }
            elsif (!exists($ambig_map{$ensembl_gene_id})) {
                push @{$ambig_map{$ensembl_gene_id}},
                    keys %{$ensembl2gene_grch37_map->{$ensembl_stable_gene_id}};
            }
        }
        $unmapped_ids{$ensembl_gene_id}++ unless exists $id_map{$ensembl_gene_id};
    }
    for my $ensembl_gene_id (natsort keys %unmapped_ids) {
        (my $ensembl_stable_gene_id = $ensembl_gene_id) =~ s/\.\d+$//;
        if (exists($ensembl2gene_ids_map->{$ensembl_stable_gene_id})) {
            my %ensembl_mapped_symbols;
            for my $entrez_gene_id (nsort keys %{$ensembl2gene_ids_map->{$ensembl_stable_gene_id}}) {
                my $mapped_symbol = find_entrez_gene_id_symbol($entrez_gene_id);
                if (defined($mapped_symbol) and !exists($mapped_ids{$mapped_symbol}) and (
                    $map_id_type eq 'hgnc_symbol' or
                    ($map_id_type eq 'custom_symbol' and exists($custom_map_ids{$mapped_symbol}))
                )) {
                    $ensembl_mapped_symbols{$mapped_symbol}++;
                }
            }
            if (%ensembl_mapped_symbols) {
                if (scalar(keys %ensembl_mapped_symbols) == 1) {
                    my ($mapped_symbol) = keys %ensembl_mapped_symbols;
                    $id_map{$ensembl_gene_id} = $mapped_symbol;
                    $mapped_ids{$mapped_symbol}++;
                    delete $ambig_map{$ensembl_gene_id};
                }
                else {
                    push @{$ambig_map{$ensembl_gene_id}}, keys %ensembl_mapped_symbols;
                }
            }
        }
        delete $unmapped_ids{$ensembl_gene_id} if exists $id_map{$ensembl_gene_id};
    }
}
# symbol
elsif ($src_id_type eq 'symbol') {
    for my $symbol (@ids) {
        if ($map_id_type eq 'custom_symbol' and exists($custom_map_ids{$symbol})) {
            $id_map{$symbol} = $symbol;
            $mapped_ids{$symbol}++;
        }
        if (!exists($id_map{$symbol}) and exists($symbol2gene_ids_map->{$symbol})) {
            if (scalar(keys %{$symbol2gene_ids_map->{$symbol}}) == 1) {
                my ($entrez_gene_id) = keys %{$symbol2gene_ids_map->{$symbol}};
                my $mapped_symbol = find_entrez_gene_id_symbol($entrez_gene_id);
                if (defined($mapped_symbol) and (
                    $map_id_type eq 'hgnc_symbol' or
                    ($map_id_type eq 'custom_symbol' and exists($custom_map_ids{$mapped_symbol}))
                )) {
                    $id_map{$symbol} = $mapped_symbol;
                    $mapped_ids{$mapped_symbol}++;
                }
            }
            else {
                for my $entrez_gene_id (nsort keys %{$symbol2gene_ids_map->{$symbol}}) {
                    if (scalar(keys %{$gene_id2symbols_map->{$entrez_gene_id}}) == 1) {
                        my $mapped_symbol = find_entrez_gene_id_symbol($entrez_gene_id);
                        if (defined($mapped_symbol) and (
                            $map_id_type eq 'hgnc_symbol' or
                            ($map_id_type eq 'custom_symbol' and exists($custom_map_ids{$mapped_symbol}))
                        )) {
                            $id_map{$symbol} = $mapped_symbol;
                            $mapped_ids{$mapped_symbol}++;
                            last;
                        }
                    }
                }
            }
        }
        if (!exists($id_map{$symbol}) and exists($gene_history_by_symbol_hashref->{$symbol})) {
            if (scalar(keys %{$gene_history_by_symbol_hashref->{$symbol}}) == 1) {
                my ($entrez_gene_id) = keys %{$gene_history_by_symbol_hashref->{$symbol}};
                my $mapped_symbol = find_entrez_gene_id_symbol($entrez_gene_id);
                if (defined($mapped_symbol) and (
                    $map_id_type eq 'hgnc_symbol' or
                    ($map_id_type eq 'custom_symbol' and exists($custom_map_ids{$mapped_symbol}))
                )) {
                    $id_map{$symbol} = $mapped_symbol;
                    $mapped_ids{$mapped_symbol}++;
                }
            }
            else {
                for my $entrez_gene_id (nsort keys %{$gene_history_by_symbol_hashref->{$symbol}}) {
                    if (scalar(keys %{$gene_id2symbols_map->{$entrez_gene_id}}) == 1) {
                        my $mapped_symbol = find_entrez_gene_id_symbol($entrez_gene_id);
                        if (defined($mapped_symbol) and (
                            $map_id_type eq 'hgnc_symbol' or
                            ($map_id_type eq 'custom_symbol' and exists($custom_map_ids{$mapped_symbol}))
                        )) {
                            $id_map{$symbol} = $mapped_symbol;
                            $mapped_ids{$mapped_symbol}++;
                            last;
                        }
                    }
                }
            }
        }
        $unmapped_ids{$symbol}++ unless exists $id_map{$symbol};
    }
    for my $symbol (natsort keys %unmapped_ids) {
        if (exists($symbol2gene_ids_map->{$symbol})) {
            my %symbol_mapped_symbols;
            for my $entrez_gene_id (nsort keys %{$symbol2gene_ids_map->{$symbol}}) {
                my $mapped_symbol = find_entrez_gene_id_symbol($entrez_gene_id);
                if (defined($mapped_symbol) and !exists($mapped_ids{$mapped_symbol}) and (
                    $map_id_type eq 'hgnc_symbol' or
                    ($map_id_type eq 'custom_symbol' and exists($custom_map_ids{$mapped_symbol}))
                )) {
                    $symbol_mapped_symbols{$mapped_symbol}++;
                }
            }
            if (%symbol_mapped_symbols) {
                if (scalar(keys %symbol_mapped_symbols) == 1) {
                    my ($mapped_symbol) = keys %symbol_mapped_symbols;
                    $id_map{$symbol} = $mapped_symbol;
                    $mapped_ids{$mapped_symbol}++;
                }
                else {
                    push @{$ambig_map{$symbol}}, keys %symbol_mapped_symbols;
                }
            }
        }
        if (!exists($id_map{$symbol}) and exists($gene_history_by_symbol_hashref->{$symbol})) {
            my %symbol_mapped_symbols;
            for my $entrez_gene_id (nsort keys %{$gene_history_by_symbol_hashref->{$symbol}}) {
                my $mapped_symbol = find_entrez_gene_id_symbol($entrez_gene_id);
                if (defined($mapped_symbol) and !exists($mapped_ids{$mapped_symbol}) and (
                    $map_id_type eq 'hgnc_symbol' or
                    ($map_id_type eq 'custom_symbol' and exists($custom_map_ids{$mapped_symbol}))
                )) {
                    $symbol_mapped_symbols{$mapped_symbol}++;
                }
            }
            if (%symbol_mapped_symbols) {
                if (scalar(keys %symbol_mapped_symbols) == 1) {
                    my ($mapped_symbol) = keys %symbol_mapped_symbols;
                    $id_map{$symbol} = $mapped_symbol;
                    $mapped_ids{$mapped_symbol}++;
                }
                else {
                    push @{$ambig_map{$symbol}}, keys %symbol_mapped_symbols;
                }
            }
        }
        delete $unmapped_ids{$symbol} if exists $id_map{$symbol};
    }
}
open(my $out_fh, '>', "${file_basename}_${out_file_suffix}${file_ext}") or die "$!";
for my $id (@ids) {
    print $out_fh "$id\t", (exists($id_map{$id}) ? $id_map{$id} : ''), "\n";
}
close($out_fh);
print scalar(keys %id_map), ' / ',
      scalar($find_entrez_map ? @ids : keys %custom_map_ids), " IDs mapped\n";
open(my $un_fh, '>', 'unmapped_ids.out') or die "$!";
print $un_fh "$_\n" for natsort keys %unmapped_ids;
close($un_fh);
if (%ambig_map) {
    for my $id (natsort keys %ambig_map) {
        print "Ambiguous map: $id => ", join(', ', natsort @{$ambig_map{$id}}), "\n";
    }
}

sub find_entrez_gene_id_symbol {
    my ($entrez_gene_id) = @_;
    $entrez_gene_id =~ s/\D//g;
    die "Invalid Gene ID: $entrez_gene_id" if $entrez_gene_id =~ /^\s*$/;
    if (exists($gene_info_hashref->{$entrez_gene_id})) {
        if ($find_entrez_map or exists($custom_map_ids{$gene_info_hashref->{$entrez_gene_id}->{symbol}})) {
            return $gene_info_hashref->{$entrez_gene_id}->{symbol};
        }
        else {
            my @custom_symbol_synonyms;
            for my $symbol_synonym (
                grep {
                    $_ ne $gene_info_hashref->{$entrez_gene_id}->{symbol}
                } keys %{$gene_id2symbols_map->{$entrez_gene_id}}
            ) {
                push @custom_symbol_synonyms, $symbol_synonym if exists($custom_map_ids{$symbol_synonym});
            }
            if (scalar(@custom_symbol_synonyms) == 1) {
                return $custom_symbol_synonyms[0];
            }
            elsif (scalar(@custom_symbol_synonyms) > 1) {
                for my $custom_symbol_synonym (@custom_symbol_synonyms) {
                    if (
                        exists($gene_symbol2id_bestmap->{$custom_symbol_synonym}) and
                        !exists($gene_symbol2id_bestmap->{$custom_symbol_synonym}->{ambig_gene_map})
                    ) {
                        return $custom_symbol_synonym;
                    }
                }
            }
            return undef;
        }
    }
    elsif (exists($gene_history_hashref->{$entrez_gene_id})) {
        if (exists($gene_history_hashref->{$entrez_gene_id}->{current_gene_id})) {
            my $current_gene_id = $gene_history_hashref->{$entrez_gene_id}->{current_gene_id};
            if ($find_entrez_map or exists($custom_map_ids{$gene_info_hashref->{$current_gene_id}->{symbol}})) {
                return $gene_info_hashref->{$current_gene_id}->{symbol};
            }
            else {
                my @custom_symbol_synonyms;
                for my $symbol_synonym (
                    grep {
                        $_ ne $gene_info_hashref->{$current_gene_id}->{symbol}
                    } keys %{$gene_id2symbols_map->{$current_gene_id}}
                ) {
                    push @custom_symbol_synonyms, $symbol_synonym if exists($custom_map_ids{$symbol_synonym});
                }
                if (scalar(@custom_symbol_synonyms) == 1) {
                    return $custom_symbol_synonyms[0];
                }
                elsif (scalar(@custom_symbol_synonyms) > 1) {
                    for my $custom_symbol_synonym (@custom_symbol_synonyms) {
                        if (
                            exists($gene_symbol2id_bestmap->{$custom_symbol_synonym}) and
                            !exists($gene_symbol2id_bestmap->{$custom_symbol_synonym}->{ambig_gene_map})
                        ) {
                            return $custom_symbol_synonym;
                        }
                    }
                }
                return undef;
            }
        }
        # discontinued Gene ID
        elsif ($find_entrez_map or exists($custom_map_ids{$gene_history_hashref->{$entrez_gene_id}->{discontinued_symbol}})) {
            # print "Discontinued Symbol: $gene_history_hashref->{$entrez_gene_id}->{discontinued_symbol}\n";
            # return $gene_history_hashref->{$entrez_gene_id}->{discontinued_symbol};
            return undef;
        }
        else {
            return undef;
        }
    }
    # invalid Gene ID
    else {
        print "Invalid Gene ID: $entrez_gene_id\n";
        return undef;
    }
}

__END__

=head1 NAME

map_gene_ids.pl - Gene ID Mapper

=head1 SYNOPSIS

 map_gene_ids.pl [options] [input data or ID file path]

 Options:
    --src-id-type=<type>            Input file ID type (required)
    --map-id-type=<type>            Output file map ID type (optional)
    --custom-symbol-file=<file>     Custom symbol map file (optional)
    --help                          Display usage and exit
    --version                       Display program version and exit

=cut
