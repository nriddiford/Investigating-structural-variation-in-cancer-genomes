#!/usr/bin/perl
use strict;
use warnings;

use feature qw /say /;

# open my $table, '>', 'markdown_table.txt' or die $!;

my $count = 0;
my @cols;
while(<DATA>){
	chomp;
	my @split = split;
	
	s/\"//g foreach @split;
	# Print header and formatting row
	if ($. == 1) {
		print "| $_ " foreach @split;
		print "|\n";
		print "|:---:" foreach @split;
		print "|\n";
	}
	
	elsif ($. > 1) {
		print "| ";
		print join(" | ", @split);
		print " |\n";
	 }

}



__DATA__
Chrom	Position	Ref	Cons	Reads1	Reads2	VarFreq	Strands1	Strands2	Qual1	Qual2	Pvalue	MapQual1	MapQual2	Reads1Plus	Reads1Minus	Reads2Plus	Reads2Minus	VarAllele
2L	5355	C	Y	21	2	8.7%	2	1	46	69	0.98	1	1	14	7	2	0	T
2L	5372	T	W	6	9	60%	2	2	60	44	0.98	1	1	5	1	8	1	A
2L	5390	T	W	13	11	45.83%	2	2	47	40	0.98	1	1	12	1	10	1	A