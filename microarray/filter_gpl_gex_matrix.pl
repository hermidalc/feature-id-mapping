#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename qw(fileparse);
use Storable qw(lock_retrieve);
use Data::Dumper;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Terse = 1;
$Data::Dumper::Deepcopy = 1;

my %id2symbol_map = %{lock_retrieve("$ARGV[0].pls")};
my ($file_basename, $file_dir, $file_ext) = fileparse($ARGV[1], qr/\.[^.]*/);
my (@col_headers);
open(my $out, '>', "${file_basename}_gene${file_ext}") or die "$!";
open(my $in, '<', $ARGV[1]) or die "$!";
while (<$in>) {
    next if m/^!/ or m/^\s*$/;
    chomp;
    if (!@col_headers and m/^"ID_REF"\t/) {
        @col_headers = split /\t/;
        s/^"|"$//g for @col_headers;
        print $out join("\t", ('', @col_headers[1..$#col_headers])), "\n";
    }
    else {
        my @fields = split /\t/;
        s/^"|"$//g for @fields;
        $fields[0] =~ s/\s+//g;
        print $out join("\t", ($id2symbol_map{$fields[0]}, @fields[1..$#fields])), "\n";
    }
}
