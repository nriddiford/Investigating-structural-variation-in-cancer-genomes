#!/usr/bin/perl
use strict;
use warnings;

use feature qw/ say /;


my $raw_cnv_counts = $ARGV[0];
my $cnv_seq_calls = $ARGV[1];

die usage() unless $raw_cnv_counts and $cnv_seq_calls;

open my $calls, '<', $cnv_seq_calls or die $!;


open my $calls_out, '>', 'cnv-seq_calls.txt' or die $!;

my $count = 0;
while(<$calls>){
	chomp;
	next if /^cnv/;
	my ($chr, $start, $stop, $log2) = (split)[1,2,3,5];
	$chr =~ s/chr//;
	print $calls_out join("\t", $chr, $start, $stop, $log2 ) . "\n";
	$count++;
}




open my $counts, '<', $raw_cnv_counts or die $!;

open my $up_scatter, '>', 'CNV_up.txt' or die $!;
open my $down_scatter, '>', 'CNV_down.txt' or die $!;

while(<$counts>){
	chomp;
	s/\"//g;
	next if /^chromosome/;
	my ($chr, $start, $stop, $log2) = (split)[0,1,2,6];
	
	next if $log2 eq 'NA';
	next if $log2 =~ /-?Inf/;
	
	my ($round) = sprintf "%.2f", $log2;
	print $up_scatter join("\t", $chr, $start, $stop, $log2 ) . "\n" if $round > 0;
	print $down_scatter join("\t", $chr, $start, $stop, $log2 ) . "\n" if $round < 0;
}

say "---------------------------------";
say " $count CNV regions called by CNV-Seq";
say "---------------------------------";

say "Writing copy number calls to 'cnv-seq_calls.txt' [CNV-Seq heatmap file]";
say "Writing copy number gains to 'CNV_up.txt' [CNV-Seq scatter file]";
say "Writing copy number losses to 'CNV_down.txt' [CNV-Seq scatter file]";


sub usage {
	say "usage $0 <.cnv> <_cnvs.txt>";
}
