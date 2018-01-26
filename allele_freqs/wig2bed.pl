#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use feature qw/ say /;

my ($chrom, $len);

my $file = shift;

open my $in, '<', $file;

while(<$in>){
  chomp;
  
  if ($_ =~ /^var/){
    ($chrom) = $_ =~ /chrom=(.*?)\s/;
    $_ =~ /span=(\d+)/ ? ($len) = $1 : $len = 1; 
    next;
  }
  next if $len < 50;
  my ($start, $score) = split;
  $start = $start - 1;
  my $end = $start + $len;
  next if not $score > 0.5;
  print join("\t", $chrom, $start, $end, $score) . "\n";
}

