#!/usr/bin/perl

use strict;
use warnings;
use feature qw/ say /;
use Data::Dumper;

use File::Basename;

exit usage() unless $#ARGV == 0;

my $file = $ARGV[0];


my ($name, $extention) = split(/\.([^.]+)$/, basename($file), 2);
my ($id) = split(/\./, $name, 2);

open my $in, '<', $file or die $!;

open my $out, '>', $id . "_control-freec" . ".gff3" or die $!;

# #Print headers to both files
print $out "##gff-version 3\n";
print $out "#track name=\"$id Freec\" gffTags=on\n";

while(<$in>){
	chomp;
	my ($chromosome, $start, $stop, $copy_no, $type) = split;
	my ($up, $down) = ("#25CAA2", "#CA254D");
	# gff if 1-based
        $start++;
        $stop++;
	if ($type eq 'gain'){
		print $out join ("\t", $chromosome, "CNV-Seq", "DUP", $start, $stop, ".", "+", ".", "Name=$type;copy_no=$copy_no;colour=$up;") . "\n";
	}
	elsif ($type eq 'loss'){
		print $out join ("\t", $chromosome, "CNV-Seq", "DEL", $start, $stop, ".", "+", ".", "Name=$type;copy_no=$copy_no;colour=$down;") . "\n";
	}
}



sub usage {
	say "Usage: $0 <_CNV file>";
}
__DATA__
2L	20150000	20200000	6	gain
2R	750000	800000	1	loss
2R	1950000	2000000	3	gain
2R	2800000	2850000	4	gain
2R	3050000	3100000	4	gain
2R	4200000	4250000	3	gain
2R	10500000	10600000	3	gain
3L	25150000	25200000	4	gain
3L	26350000	26400000	4	gain
3L	26950000	27000000	1	loss
3L	27250000	27300000	3	gain
3L	27750000	27800000	3	gain
3R	1350000	1400000	5	gain
3R	2050000	2100000	10	gain
3R	2350000	2400000	3	gain
X	2800000	3250000	0	loss
X	14700000	14750000	2	gain
X	19650000	19700000	2	gain
X	22900000	23542271	2	gain
X	22900000	23542271	2	gain

