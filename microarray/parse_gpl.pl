#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename qw(fileparse);
use Storable qw(lock_nstore);
use Data::Dumper;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Terse = 1;
$Data::Dumper::Deepcopy = 1;

my $gpl_file = $ARGV[0];
my (@col_headers, %col_header_idxs, %id2symbol_map);
open(my $in, '<', $gpl_file) or die "$!";
while (<$in>) {
    next if m/^#/;
    chomp;
    if (!@col_headers and m/^ID\t/) {
        @col_headers = split /\t/;
        %col_header_idxs = map { $col_headers[$_] => $_ } 0 .. $#col_headers;
    }
    else {
        my @fields = split /\t/;
        my $symbol_idx = exists($col_header_idxs{'Gene Symbol'})
                       ? $col_header_idxs{'Gene Symbol'}
                       : exists($col_header_idxs{'Symbol'})
                           ? $col_header_idxs{'Symbol'}
                           : exists($col_header_idxs{'gene_assignment'})
                               ? $col_header_idxs{'gene_assignment'}
                               : die "Doh\n";
        my $symbol = exists($col_header_idxs{'gene_assignment'})
                   ? $fields[$symbol_idx] ne "---"
                       ? @{[split /\/\//, $fields[$symbol_idx]]}[1]
                       : ''
                   : $fields[$symbol_idx] =~ m/\/\/\//
                       ? @{[split /\/\/\//, $fields[$symbol_idx]]}[0]
                       :  $fields[$symbol_idx];
        $symbol =~ s/\s+//g;
        $symbol = '' if $symbol eq '---';
        $id2symbol_map{$fields[$col_header_idxs{'ID'}]} = $symbol;
    }
}
close($in);
my ($file_basename, $file_dir, $file_ext) = fileparse($gpl_file, qr/\.[^.]*/);
lock_nstore(\%id2symbol_map, "$file_basename.pls");
print Dumper(\%id2symbol_map);
