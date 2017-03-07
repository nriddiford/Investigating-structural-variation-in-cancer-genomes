#!/usr/bin/perl

use strict;
use warnings;

use feature qw/ say /;

my @files = qw/ lumpy_svs.txt lumpy_bnds.txt cnv-seq_calls.txt delly_bnds.txt delly_svs.txt freec_calls.txt /;

for (@files){
	open my $file, '<', $_ or die $!;
	say "$_";
}