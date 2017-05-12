#/usr/bin/perl
use strict;
use warnings;
use feature qw /say/;

my @files = @ARGV;

my @keys = qw / 2L 2R 3L 3R 4 X Y /;
my %filter;
$filter{$_} = 1 for (@keys);

for my $f (@files){
open my $in, '<', $f or die $!;
open my $out, '>', "$f\.filt" or die $!;

    while(<$in>){
        chomp;
        my @split = split;
        print $out "$_\n" if $filter{$split[0]};
    }
}
