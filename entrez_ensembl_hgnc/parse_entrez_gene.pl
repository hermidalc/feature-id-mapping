#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use File::Basename qw(fileparse);
use File::Fetch;
use Getopt::Long qw(:config auto_help auto_version);
use Hash::Util;
use Pod::Usage qw(pod2usage);
use Sort::Key qw(nsort);
use Sort::Key::Natural qw(natsort rnatsort natkeysort);
use Storable qw(lock_nstore);
use Data::Dumper;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Terse = 1;
$Data::Dumper::Deepcopy = 1;

# unbuffer error and output streams (make sure STDOUT is last so that it remains the default filehandle)
select(STDERR); $| = 1;
select(STDOUT); $| = 1;

sub is_integer {
    my $value = shift;
    return $value =~ m/^-?\d+$/o ? 1 : 0;
}

# config
my $NCBI_FTP_BASE_URI = 'ftp://ftp.ncbi.nih.gov';
my $ENTREZ_GENE_DATA_BASE_URI = "$NCBI_FTP_BASE_URI/gene/DATA";
my %CTK_ENTREZ_GENE_DATA = (
    # 'gene_info' => {
    #     file_uri => "$ENTREZ_GENE_DATA_BASE_URI/gene_info.gz",
    # },
    'gene_history' => {
        file_uri => "$ENTREZ_GENE_DATA_BASE_URI/gene_history.gz",
    },
    # 'gene2refseq' => {
    #     file_uri => "$ENTREZ_GENE_DATA_BASE_URI/gene2refseq.gz",
    # },
    # 'gene2accession' => {
    #     file_uri => "$ENTREZ_GENE_DATA_BASE_URI/gene2accession.gz",
    # },
    'gene2ensembl' => {
        file_uri => "$ENTREZ_GENE_DATA_BASE_URI/gene2ensembl.gz",
    },
    # 'gene2unigene' => {
    #     # NCBI doesn't gzip gene2unigene because file is small
    #     file_uri => "$ENTREZ_GENE_DATA_BASE_URI/gene2unigene",
    # },
);
my %CTK_ENTREZ_GENE_ORGANISM_DATA = (
    'Homo sapiens' => {
        gene_info_file_uri => "$ENTREZ_GENE_DATA_BASE_URI/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
        tax_id => '9606',
    },
    # 'Mus musculus' => {
    #     gene_info_file_uri => "$ENTREZ_GENE_DATA_BASE_URI/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz",
    #     tax_id => '10090',
    # },
    # 'Rattus norvegicus' => {
    #     gene_info_file_uri => "$ENTREZ_GENE_DATA_BASE_URI/GENE_INFO/Mammalia/Rattus_norvegicus.gene_info.gz",
    #     tax_id => '10116',
    # },
);
my $CTK_DATA_ID_MAPPING_GENE_SYMBOL_SUFFIX = '_gene_symbols';
# --entrez-download download latest Entrez Gene data (default false)
my $entrez_download = 0;
my $debug = 0;
my $verbose = 0;
GetOptions(
    'entrez-download' => \$entrez_download,
    'debug'           => \$debug,
    'verbose'         => \$verbose,
) || pod2usage(-verbose => 0);
# entrez gene
my $entrez_data_work_dir = $FindBin::Bin;
my $entrez_data_tmp_dir = $FindBin::Bin;
my (%all_tax_ids, $gene_info_hashref, $gene_history_hashref, $gene_history_by_symbol_hashref,
    $symbol2gene_ids_map, $uc_symbol2gene_ids_map, $gene_id2symbols_map, $gene_symbol2id_bestmap,
    $accession2gene_ids_map, $ensembl2gene_ids_map, $gene2ensembl_ids_map, $unigene2gene_ids_map);
print "[Entrez Gene Data]\n";
if ($entrez_download) {
    # download and uncompress latest public master Entrez Gene data files
    $entrez_data_work_dir = $entrez_data_tmp_dir;
    for my $entrez_gene_file_key (rnatsort keys %CTK_ENTREZ_GENE_DATA) {
        my $file_uri = $CTK_ENTREZ_GENE_DATA{$entrez_gene_file_key}{file_uri};
        my $ff = File::Fetch->new(uri => $file_uri) or die "\n\nERROR: File::Fetch object constructor error\n\n";
        my ($gi_file_basename, undef, $gi_file_ext) = fileparse($file_uri, qr/\.[^.]*/);
        print "Fetching latest $gi_file_basename$gi_file_ext\n";
        $ff->fetch(to => $entrez_data_work_dir) or die "\n\nERROR: File::Fetch fetch error: ", $ff->error, "\n\n";
        print "Uncompressing $gi_file_basename$gi_file_ext\n";
        if ($gi_file_ext) {
            my $uncompress_cmd = lc($gi_file_ext) eq '.gz'  ? "gzip -df $entrez_data_work_dir/$gi_file_basename$gi_file_ext"
                               : lc($gi_file_ext) eq '.zip' ? "unzip -oq $entrez_data_work_dir/$gi_file_basename$gi_file_ext -d $entrez_data_work_dir"
                               : die "\n\nERROR: unsupported compressed file extension '$gi_file_ext'\n\n";
            system($uncompress_cmd) == 0 or die "\nERROR: $uncompress_cmd system call error: ", $? >> 8, "\n\n";
        }
    }
    # download and uncompress organism-specific gene_info files or parse master gene_info file to create organism-specific if one doesn't exist at NCBI
    for my $organism_name (sort keys %CTK_ENTREZ_GENE_ORGANISM_DATA) {
        if (exists $CTK_ENTREZ_GENE_ORGANISM_DATA{$organism_name}{gene_info_file_uri} and
            defined $CTK_ENTREZ_GENE_ORGANISM_DATA{$organism_name}{gene_info_file_uri}) {
            my $ff = File::Fetch->new(uri => $CTK_ENTREZ_GENE_ORGANISM_DATA{$organism_name}{gene_info_file_uri})
                or die "\n\nERROR: File::Fetch object constructor error\n\n";
            my ($gi_file_basename, undef, $gi_file_ext) = fileparse($CTK_ENTREZ_GENE_ORGANISM_DATA{$organism_name}{gene_info_file_uri}, qr/\.[^.]*/);
            print "Fetching latest $gi_file_basename$gi_file_ext\n";
            $ff->fetch(to => $entrez_data_work_dir) or die "\n\nERROR: File::Fetch fetch error: ", $ff->error, "\n\n";
            print "Uncompressing $gi_file_basename$gi_file_ext\n";
            if ($gi_file_ext) {
                my $uncompress_cmd = lc($gi_file_ext) eq '.gz'  ? "gzip -df $entrez_data_work_dir/$gi_file_basename$gi_file_ext"
                                   : lc($gi_file_ext) eq '.zip' ? "unzip -oq $entrez_data_work_dir/$gi_file_basename$gi_file_ext -d $entrez_data_work_dir"
                                   : die "\n\nERROR: unsupported compressed file extension '$gi_file_ext'\n\n";
                system($uncompress_cmd) == 0 or die "\nERROR: $uncompress_cmd system call error: ", $? >> 8, "\n\n";
            }
        }
        elsif (exists $CTK_ENTREZ_GENE_ORGANISM_DATA{$organism_name}{gene_info_tax_ids} and
               defined $CTK_ENTREZ_GENE_ORGANISM_DATA{$organism_name}{gene_info_tax_ids}) {
            (my $organism_file_basename = $organism_name) =~ s/\s+/_/g;
            print "Parsing master gene_info file to create $organism_file_basename.gene_info, ",
                  "using taxonomy IDs ", join(',', @{$CTK_ENTREZ_GENE_ORGANISM_DATA{$organism_name}{gene_info_tax_ids}}), ': ';
            my %tax_ids = map { $_ => 1 } @{$CTK_ENTREZ_GENE_ORGANISM_DATA{$organism_name}{gene_info_tax_ids}};
            open(my $gene_info_fh, '<', "$entrez_data_work_dir/gene_info")
                or die "\n\nERROR: could not open $entrez_data_work_dir/gene_info: $!\n\n";
            open(my $output_fh, '>', "$entrez_data_work_dir/$organism_file_basename.gene_info")
                or die "\n\nERROR: could not create $entrez_data_work_dir/$organism_file_basename.gene_info: $!\n\n";
            # skip column header for organism-specific gene_info file
            my $col_header = <$gene_info_fh>;
            my $genes_added = 0;
            while (<$gene_info_fh>) {
                my ($tax_id) = split /\t/;
                if ($tax_ids{$tax_id}) {
                    print $output_fh $_;
                    $genes_added++;
                }
            }
            close($output_fh);
            close($gene_info_fh);
            print "$genes_added genes added\n";
            die "\n\nERROR: problem with parsing or tax IDs, 0 genes added!\n\n" if $genes_added == 0;
        }
        else {
            die "\n\nERROR: problem with $organism_name Entrez Gene configuration! ",
                "Missing gene_info_file_uri or gene_info_tax_ids, please check the configuration file.\n\n";
        }
    }
}
# parse latest Entrez Gene data from organism-specific gene_info files, create data structures and gene symbol mappings
my $total_genes_parsed = 0;
for my $organism_name (sort keys %CTK_ENTREZ_GENE_ORGANISM_DATA) {
    my %organism_gene_symbols;
    my $genes_parsed = 0;
    (my $organism_file_basename = $organism_name) =~ s/\s+/_/g;
    print "Parsing $organism_file_basename.gene_info: ";
    open(my $gene_info_fh, '<', "$entrez_data_work_dir/$organism_file_basename.gene_info")
        or die "ERROR: could not open $entrez_data_work_dir/$organism_file_basename.gene_info: $!";
    while (<$gene_info_fh>) {
        m/^#/ && next;
        my ($tax_id, $gene_id, $gene_symbol, undef, $symbol_synonyms_str, undef, $chromosome, undef, $gene_desc) = split /\t/;
        s/\s+//g for $tax_id, $gene_id, $gene_symbol, $chromosome;
        $symbol_synonyms_str =~ s/^\s+//;
        $symbol_synonyms_str =~ s/\s+$//;
        die "ERROR: Organism Tax ID $tax_id is not an integer" unless is_integer($tax_id);
        die "ERROR: Gene ID $gene_id is not an integer" unless is_integer($gene_id);
        die "ERROR: Gene ID $gene_id defined more than once (there is a problem with Entrez Gene gene_info file)" if exists $gene_info_hashref->{$gene_id};
        $all_tax_ids{$tax_id}++;
        $gene_info_hashref->{$gene_id}->{organism_tax_id} = $tax_id;
        $gene_info_hashref->{$gene_id}->{symbol} = $gene_symbol;
        $gene_info_hashref->{$gene_id}->{description} = ($gene_desc and $gene_desc ne '-') ? $gene_desc : undef;
        $gene_info_hashref->{$gene_id}->{synonyms} = ($symbol_synonyms_str and $symbol_synonyms_str ne '-') ? $symbol_synonyms_str : undef;
        # add gene symbols and maps for organism and only those from organism mitochondria which don't already exist for organism
        if ((exists $CTK_ENTREZ_GENE_ORGANISM_DATA{$organism_name}{tax_id} and
             defined $CTK_ENTREZ_GENE_ORGANISM_DATA{$organism_name}{tax_id} and
             $tax_id == $CTK_ENTREZ_GENE_ORGANISM_DATA{$organism_name}{tax_id}) or
            !exists $organism_gene_symbols{$gene_symbol}) {
            $organism_gene_symbols{$gene_symbol}++;
            $symbol2gene_ids_map->{$organism_name}->{$gene_symbol}->{$gene_id}++;
            $uc_symbol2gene_ids_map->{$organism_name}->{uc($gene_symbol)}->{$gene_id}++;
            $gene_id2symbols_map->{$organism_name}->{$gene_id}->{$gene_symbol}++;
        }
        $genes_parsed++;
    }
    close($gene_info_fh);
    # organism gene symbol synonyms need to be processed after loading all organism official gene symbols because
    # I need to know all official gene symbols to determine which gene symbol synonyms are really valid or not
    for my $gene_id (keys %{$gene_id2symbols_map->{$organism_name}}) {
        next unless defined $gene_info_hashref->{$gene_id}->{synonyms};
        my @valid_symbol_synonyms;
        for my $symbol_synonym (split /\|/, $gene_info_hashref->{$gene_id}->{synonyms}) {
            $symbol_synonym =~ s/^\s+//;
            $symbol_synonym =~ s/\s+$//;
            # if gene symbol synonym is not the same as an official gene symbol for this organism then it is valid, otherwise invalid
            if (!exists $organism_gene_symbols{$symbol_synonym}) {
                $symbol2gene_ids_map->{$organism_name}->{$symbol_synonym}->{$gene_id}++;
                $uc_symbol2gene_ids_map->{$organism_name}->{uc($symbol_synonym)}->{$gene_id}++;
                # reverse map of gene symbol synonyms doesn't really make sense but it's only for reference
                $gene_id2symbols_map->{$organism_name}->{$gene_id}->{$symbol_synonym}++;
                # to update $gene_info_hashref synonyms below
                push @valid_symbol_synonyms, $symbol_synonym;
            }
        }
        $gene_info_hashref->{$gene_id}->{synonyms} = join '|', @valid_symbol_synonyms;
    }
    print "$genes_parsed genes\n";
    die "ERROR: problem parsing, 0 genes parsed!" if $genes_parsed == 0;
    $total_genes_parsed += $genes_parsed;
}
# print "--> $total_genes_parsed <-- total current genes\n";
# parse latest Entrez Gene gene_history
my $relevant_old_genes_processed = 0;
print "Parsing gene_history: ";
open(my $gene_history_fh, '<', "$entrez_data_work_dir/gene_history") or die "ERROR: Could not open $entrez_data_work_dir/gene_history: $!";
while (<$gene_history_fh>) {
    m/^#/ && next;
    my ($tax_id, $current_gene_id, $discontinued_gene_id, $discontinued_symbol) = split /\t/;
    s/\s+//g for $tax_id, $current_gene_id, $discontinued_gene_id, $discontinued_symbol;
    die "ERROR: Organism Tax ID $tax_id is not an integer" unless is_integer($tax_id);
    die "ERROR: Gene ID $current_gene_id is not an integer" unless $current_gene_id eq '-' or is_integer($current_gene_id);
    die "ERROR: discontinued Gene ID $discontinued_gene_id is not an integer" unless is_integer($discontinued_gene_id);
    die "ERROR: discontinued Gene ID $discontinued_gene_id defined in current Entrez Gene data (there is a problem with Entrez Gene)"
        if exists $gene_info_hashref->{$discontinued_gene_id};
    # only load gene history data lines for organisms used in our installation
    if (exists $all_tax_ids{$tax_id}) {
        if (!exists $gene_history_hashref->{$discontinued_gene_id}) {
            $gene_history_hashref->{$discontinued_gene_id}->{organism_tax_id} = $tax_id;
            if ($current_gene_id and $current_gene_id ne '-') {
                $gene_history_hashref->{$discontinued_gene_id}->{current_gene_id} = $current_gene_id;
            }
            if ($discontinued_symbol and $discontinued_symbol ne '-') {
                $gene_history_hashref->{$discontinued_gene_id}->{discontinued_symbol} = $discontinued_symbol;
            }
            if (
                $current_gene_id and $current_gene_id ne '-' and
                $discontinued_symbol and $discontinued_symbol ne '-'
            ) {
                $gene_history_by_symbol_hashref->{$discontinued_symbol}->{$current_gene_id}++;
            }
        }
        else {
            print "ERROR: discontinued Gene ID $discontinued_gene_id defined more than once (there is a problem with Entrez Gene)\n";
        }
        $relevant_old_genes_processed++;
    }
}
close($gene_history_fh);
print "$relevant_old_genes_processed relevant discontinued gene IDs processed\n";
# parse latest Entrez Gene gene2ensembl
my $ensembl_ids_processed = 0;
print "Parsing gene2ensembl: ";
open(my $gene2ensembl_fh, '<', "$entrez_data_work_dir/gene2ensembl") or die "ERROR: Could not open $entrez_data_work_dir/gene2ensembl: $!";
while (<$gene2ensembl_fh>) {
    m/^#/ && next;
    my ($tax_id, $gene_id, $ensembl_gene_id) = split /\t/;
    s/\s+//g for $tax_id, $gene_id, $ensembl_gene_id;
    # skip data lines for organisms not used in our installation
    next unless exists $gene_info_hashref->{$gene_id};
    die "\n\nERROR: Organism Tax ID $tax_id is not an integer\n\n" unless is_integer($tax_id);
    die "\n\nERROR: Gene ID $gene_id is not an integer\n\n" unless is_integer($gene_id);
    $ensembl2gene_ids_map->{$ensembl_gene_id}->{$gene_id}++;
    $gene2ensembl_ids_map->{$gene_id}->{$ensembl_gene_id}++;
    $ensembl_ids_processed++;
}
close($gene2ensembl_fh);
print "$ensembl_ids_processed relevant ensembl entries processed\n";
# Hash::Util::lock_hashref_recurse is broken in Perls < 5.14
if ($^V >= 'v5.14') {
    # lock data structures
    Hash::Util::lock_hashref_recurse($gene_info_hashref);
    Hash::Util::lock_hashref_recurse($gene_history_hashref);
    Hash::Util::lock_hashref_recurse($gene_history_by_symbol_hashref);
    Hash::Util::lock_hashref_recurse($symbol2gene_ids_map);
    Hash::Util::lock_hashref_recurse($uc_symbol2gene_ids_map);
    Hash::Util::lock_hashref_recurse($gene_id2symbols_map);
    Hash::Util::lock_hashref_recurse($accession2gene_ids_map);
    Hash::Util::lock_hashref_recurse($ensembl2gene_ids_map);
    Hash::Util::lock_hashref_recurse($gene2ensembl_ids_map);
    Hash::Util::lock_hashref_recurse($unigene2gene_ids_map);
}
# serialize and store data structures
print "Serializing and storing gene_info.pls\n";
lock_nstore($gene_info_hashref, "$entrez_data_tmp_dir/gene_info.pls")
    or die "ERROR: could not serialize and store to $entrez_data_tmp_dir/gene_info.pls: $!";
print "Serializing and storing gene_history.pls\n";
lock_nstore($gene_history_hashref, "$entrez_data_tmp_dir/gene_history.pls")
    or die "ERROR: could not serialize and store to $entrez_data_tmp_dir/gene_history.pls: $!";
print "Serializing and storing gene_history_by_symbol.pls\n";
lock_nstore($gene_history_by_symbol_hashref, "$entrez_data_tmp_dir/gene_history_by_symbol.pls")
    or die "ERROR: could not serialize and store to $entrez_data_tmp_dir/gene_history_by_symbol.pls: $!";
print "Serializing and storing ensembl2gene_ids.pls\n";
lock_nstore($ensembl2gene_ids_map, "$entrez_data_tmp_dir/ensembl2gene_ids.pls")
    or die "ERROR: could not serialize and store to $entrez_data_tmp_dir/ensembl2gene_ids.pls: $!";
print "Serializing and storing gene2ensembl_ids.pls\n";
lock_nstore($gene2ensembl_ids_map, "$entrez_data_tmp_dir/gene2ensembl_ids.pls")
    or die "ERROR: could not serialize and store to $entrez_data_tmp_dir/gene2ensembl_ids.pls: $!";
# generate organism gene symbol mapping files
my $mapping_data_work_dir = $FindBin::Bin;
print "[Gene Symbol Mapping]\n";
for my $organism_name (sort keys %CTK_ENTREZ_GENE_ORGANISM_DATA) {
    (my $organism_file_basename = $organism_name) =~ s/\s+/_/g;
    my $map_file_basename = "${organism_file_basename}${CTK_DATA_ID_MAPPING_GENE_SYMBOL_SUFFIX}";
    my $num_maps_written = 0;
    print "Generating ${map_file_basename}.map: ";
    open(my $map_fh, '>', "$mapping_data_work_dir/${map_file_basename}.map")
        or die "Could not create mapping file $mapping_data_work_dir/${map_file_basename}.map: $!";
    print $map_fh "Gene Symbol\tEntrez Gene IDs\n";
    for my $gene_symbol (natkeysort { lc } keys %{$symbol2gene_ids_map->{$organism_name}}) {
        my @gene_ids = nsort keys %{$symbol2gene_ids_map->{$organism_name}->{$gene_symbol}};
        print $map_fh "$gene_symbol\t", join("\t", @gene_ids), "\n";
        if (scalar(@gene_ids) == 1) {
            $gene_symbol2id_bestmap->{$organism_name}->{$gene_symbol}->{gene_id} = $gene_ids[0];
        }
        else {
            $gene_symbol2id_bestmap->{$organism_name}->{$gene_symbol}->{gene_id} = undef;
            $gene_symbol2id_bestmap->{$organism_name}->{$gene_symbol}->{ambig_gene_map}++;
        }
        $num_maps_written++;
    }
    close($map_fh);
    print "$num_maps_written maps\n";
    $num_maps_written = 0;
    print "Generating ${map_file_basename}.ucmap: ";
    open(my $ucmap_fh, '>', "$mapping_data_work_dir/${map_file_basename}.ucmap")
        or die "Could not create mapping file $mapping_data_work_dir/${map_file_basename}.ucmap: $!";
    print $ucmap_fh "Gene Symbol\tEntrez Gene IDs\n";
    for my $uc_gene_symbol (natsort keys %{$uc_symbol2gene_ids_map->{$organism_name}}) {
        my @gene_ids = nsort keys %{$uc_symbol2gene_ids_map->{$organism_name}->{$uc_gene_symbol}};
        print $ucmap_fh "$uc_gene_symbol\t", join("\t", @gene_ids), "\n";
        $num_maps_written++;
    }
    close($ucmap_fh);
    print "$num_maps_written maps\n";
    # write out reverse map file (not used by Confero only for human reference)
    $num_maps_written = 0;
    print "Generating ${map_file_basename}.revmap: ";
    open(my $revmap_fh, '>', "$mapping_data_work_dir/${map_file_basename}.revmap")
        or die "Could not create mapping file $mapping_data_work_dir/${map_file_basename}.revmap: $!";
    print $revmap_fh "Entrez Gene ID\tGene Symbols (first in list is official symbol)\n";
    for my $gene_id (nsort keys %{$gene_id2symbols_map->{$organism_name}}) {
        # first gene symbol is official symbol, followed by symbol synonyms
        print $revmap_fh "$gene_id\t", join("\t",
            $gene_info_hashref->{$gene_id}->{symbol},
            natkeysort { lc } grep { $_ ne $gene_info_hashref->{$gene_id}->{symbol} } keys %{$gene_id2symbols_map->{$organism_name}->{$gene_id}}
        ), "\n";
        $num_maps_written++;
    }
    close($revmap_fh);
    print "$num_maps_written maps\n";
    # write out best map file
    $num_maps_written = 0;
    print "Generating ${map_file_basename}.bestmap: ";
    open(my $bestmap_fh, '>', "$mapping_data_work_dir/${map_file_basename}.bestmap")
        or die "Could not create mapping file $mapping_data_work_dir/${map_file_basename}.bestmap: $!";
    print $bestmap_fh "Gene Symbol\tEntrez Gene ID\n";
    for my $gene_symbol (natkeysort { lc } keys %{$gene_symbol2id_bestmap->{$organism_name}}) {
        print $bestmap_fh "$gene_symbol\t", defined $gene_symbol2id_bestmap->{$organism_name}->{$gene_symbol}->{gene_id}
            ? $gene_symbol2id_bestmap->{$organism_name}->{$gene_symbol}->{gene_id}
            : '',
            "\n";
        $num_maps_written++;
    }
    close($bestmap_fh);
    print "$num_maps_written maps\n";
    # Hash::Util::lock_hashref_recurse is broken in Perls < 5.14
    if ($^V >= 'v5.14') {
        # lock data structure
        Hash::Util::lock_hashref_recurse($gene_symbol2id_bestmap->{$organism_name});
    }
    # serialize and store data structures
    print "Serializing and storing ${map_file_basename}.map.pls\n";
    lock_nstore($symbol2gene_ids_map->{$organism_name}, "$mapping_data_work_dir/${map_file_basename}.map.pls")
       or die "ERROR: could not serialize and store to $mapping_data_work_dir/${map_file_basename}.map.pls: $!";
    print "Serializing and storing ${map_file_basename}.ucmap.pls\n";
    lock_nstore($uc_symbol2gene_ids_map->{$organism_name}, "$mapping_data_work_dir/${map_file_basename}.ucmap.pls")
        or die "ERROR: could not serialize and store to $mapping_data_work_dir/${map_file_basename}.ucmap.pls: $!";
    print "Serializing and storing ${map_file_basename}.revmap.pls\n";
    lock_nstore($gene_id2symbols_map->{$organism_name}, "$mapping_data_work_dir/${map_file_basename}.revmap.pls")
       or die "ERROR: could not serialize and store to $mapping_data_work_dir/${map_file_basename}.revmap.pls: $!";
    print "Serializing and storing ${map_file_basename}.bestmap.pls\n";
    lock_nstore($gene_symbol2id_bestmap->{$organism_name}, "$mapping_data_work_dir/${map_file_basename}.bestmap.pls")
        or die "ERROR: could not serialize and store to $mapping_data_work_dir/${map_file_basename}.bestmap.pls: $!";
}
