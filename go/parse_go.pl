#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename qw( fileparse );
use List::Util qw( none );
use Storable qw( lock_nstore );
use Data::Dumper;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Terse = 1;
$Data::Dumper::Deepcopy = 1;

my $go_file = $ARGV[0];
my $go_symbol_data = {};
open(my $in, '<', $go_file) or die "$!";
while (<$in>) {
    chomp;
    my ($symbol, $synonyms_str, undef, $desc) = split(/\t/);
    my @synonyms = split(/\|/, $synonyms_str);
    if (
        !exists($go_symbol_data->{official}->{$symbol}) and
        none { m/^UniProtKB:/ } @synonyms
    ) {
        $go_symbol_data->{official}->{$symbol}++;
        for my $synonym (@synonyms) {
            if (!exists($go_symbol_data->{synonym}->{$synonym})) {
                $go_symbol_data->{synonym}->{$synonym} = $symbol;
            }
            # ambiguous
            else {
                # print "$synonym\n";
            }
        }
    }
    # duplicate
    else {
        # print "$_\n";
    }
}
close($in);
my ($file_basename, $file_dir, $file_ext) = fileparse($go_file, qr/\.[^.]*/);
lock_nstore($go_symbol_data, "$file_basename.pls");
print Dumper($go_symbol_data);
