#!/usr/bin/perl

use strict;
use warnings;
use feature qw/ say /;

use File::Basename;

exit usage() unless $#ARGV == 0;

my $file = $ARGV[0];

my ($name, $extention) = split(/\.([^.]+)$/, basename($file), 2);

my ($id) = split(/\_/, $name, 2);

open my $in, '<', $file or die $!;
open my $out, '>', $id . "_cnv-seq" . ".gff3" or die $!;

# Print headers to file
print $out "##gff-version 3\n";
print $out "#track name=\"$id CNV-Seq\" gffTags=on\n";

while(<$in>){
	chomp;
	next unless /^CNVR/;
	next if (split)[5] =~ /-?Inf/;

	my ($cnv, $chromosome, $start, $stop, $size, $log2, $p_value) = split;

	s/chr//g;
	my ($up, $down) = ("#25CAA2", "#CA254D");

	# gff is 1-based
        $start++;
        $stop++;

	if ($log2 > 0){
		print $out join ("\t", $chromosome, "CNV-Seq", "DUP", $start, $stop, ".", "+", ".", "Name=$cnv;log2=$log2;colour=$up;") . "\n";
	}
	elsif ($log2 < 0){
		print $out join ("\t", $chromosome, "CNV-Seq", "DEL", $start, $stop, ".", "+", ".", "Name=$cnv;log2=$log2;colour=$down;") . "\n";
	}	 
}



sub usage {
        say $0;
        say "Program to read CNVs predicted by CNV-Seq function:";
        say "'> cnv.print(data)'";
        say "Expects data in the following format:";
	say "cnv	chromosome	start	end	size	log2	p.value";
        say "Outputs a .gff3 file compatible with IGV viewing";
	say "Usage: $0 <_cnvs.txt>";
        say "Nick Riddiford 2017";
}


