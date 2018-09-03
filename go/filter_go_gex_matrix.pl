#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename qw( fileparse );
use List::Util qw( any );
use Sort::Key::Natural qw( natsort );
use Storable qw( lock_retrieve );
use Data::Dumper;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Terse = 1;
$Data::Dumper::Deepcopy = 1;

my $go_symbol_data = lock_retrieve("$ARGV[0].pls");
for my $file (@ARGV[1..$#ARGV]) {
    my ($file_basename, $file_dir, $file_ext) = fileparse($file, qr/\.[^.]*/);
    print "Parsing ${file_basename}${file_ext}\n";
    open(my $out, '>', "${file_basename}_filtered${file_ext}") or die "$!";
    open(my $in, '<', $file) or die "$!";
    my $read_header;
    my %go_symbol_found = map { $_ => 0 } keys %{$go_symbol_data->{official}};
    while (<$in>) {
        next if m/^!/ or m/^\s*$/;
        chomp;
        if (!$read_header and m/^ID_REF\t/) {
            $read_header++;
            print $out "$_\n";
        }
        else {
            my @fields = split /\t/;
            $fields[0] =~ s/\s+//g;
            if (exists($go_symbol_data->{official}->{$fields[0]})) {
                print $out "$_\n";
                $go_symbol_found{$fields[0]}++;
            }
            # elsif (exists($go_symbol_data->{synonym}->{$fields[0]})) {
            #     # print "$fields[0]"; <STDIN>;
            #     $fields[0] = $go_symbol_data->{synonym}->{$fields[0]};
            #     print $out join("\t", @fields), "\n";
            #     $go_symbol_found{$fields[0]}++;
            # }
        }
    }
    for my $symbol (natsort keys %go_symbol_found) {
        next if $go_symbol_found{$symbol};
        # print "$symbol\n";
    }
    close($out);
    close($in);
}
