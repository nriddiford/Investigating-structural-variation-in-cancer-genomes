#!/usr/bin/perl

use strict;
use warnings;

use feature qw/ say /;

my $CF_cnvs = $ARGV[0];

exit usage() unless $CF_cnvs;

open my $in, '<', $CF_cnvs or die $!;

open my $freec_out, '>', 'freec_calls.txt' or die $!;

my $cnv_count = 0;

while(<$in>){
	chomp;
	my ($chr, $start, $stop, $ploidy, $type) = (split);
	print $freec_out join("\t", $chr, $start, $stop, $ploidy ) . "\n";
	my $length = $stop - $start;
	$cnv_count++;
}

say "---------------------------------";
say " $cnv_count CNV regions called by Freec";
say "---------------------------------";

say "Writing copy number calls to 'freec_calls.txt' [Freec heatmap file]";




sub usage {
	say "Usage: $0 <_CNVS>";
}


__DATA__
2R	1685000	1694999	4	gain
2R	1690000	1704999	6	gain
2R	2265000	2274999	4	gain
2R	2270000	2284999	6	gain
2R	4310000	4324999	4	gain
3L	0	89999	2	loss
3L	5860000	5874999	5	gain
3L	24595000	24654999	4	gain
3L	25850000	25869999	4	gain
3L	26875000	27014999	4	gain
3L	27065000	27084999	4	gain
3L	27725000	27734999	4	gain
3L	27730000	27744999	6	gain
3R	1035000	1044999	5	gain
3R	1040000	1049999	7	gain
3R	1045000	1054999	4	gain
3R	4880000	4889999	5	gain
3R	4885000	4894999	9	gain
3R	4890000	4899999	7	gain
3R	8900000	8909999	4	gain
3R	8905000	8914999	7	gain
3R	8910000	8919999	5	gain
3R	31040000	31074999	4	gain
4	1200000	1209999	2	loss
X	0	1129999	1	loss
X	1125000	9949999	2	gain
X	9945000	12469999	1	loss
X	12465000	12579999	0	loss
X	12575000	15034999	2	gain
X	15030000	15044999	5	gain
X	15040000	23542271	2	gain
Y	0	269999	2	gain
Y	265000	474999	1	loss
Y	470000	479999	2	gain
Y	475000	489999	4	gain
Y	485000	1104999	2	gain
Y	1100000	1419999	1	loss
Y	1415000	1459999	2	gain
Y	1455000	3304999	1	loss
Y	3300000	3667352	2	gain
