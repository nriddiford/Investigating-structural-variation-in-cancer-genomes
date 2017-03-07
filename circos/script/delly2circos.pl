#!/usr/bin/perl

use strict;
use warnings;

use feature qw/ say /;

use Data::Dumper;

use Getopt::Long qw/ GetOptions /;

my $defined_svs;
my $vcf;
my $filter;

GetOptions( 'vcf_in=s'					=>			\$vcf,
			'defined_variants=s'		=>			\$defined_svs,
			'filter_set=s'				=>			\$filter
	  	  )	or die usage();
		  
exit usage() unless $vcf;

my %sv_id;
my $characterised_varaints;

if ($defined_svs) {
	open $characterised_varaints , '<', $defined_svs or die $!;
	while(<$characterised_varaints>){
		chomp;
		my ($id, $type) = split;
		$sv_id{$id} = $type;
	}

	say "\nAnnotated:";
	print Dumper \%sv_id;

	
}

my $filter_set;
my %id_filter;

if ($filter) {
	open $filter_set , '<', $filter or die $!;
	while(<$filter_set>){
		chomp;
		my ($id, $reason) = split;
		$id_filter{$id} = $reason;
	}
	say "\nFiltered:";
	print Dumper \%id_filter;
}

open my $in, '<', $vcf or die $!;

open my $trans_links, '>', 'delly_translocation_links.txt' or die $!;
open my $bnds, '>', 'delly_bnds.txt' or die $!;

open my $delly_svs, '>', 'delly_svs.txt' or die $!;

my @keys = qw / 2L 2R 3L 3R 4 X Y /;
my %filter;



$filter{$_} = 1 for (@keys);

my $sv_count = 0;
my $bp_count = 0;

while(<$in>){
	chomp;
	next if /^#/;
	my ($c1, $pos1, $filter) = (split)[0,1,6];
	next unless $filter{$c1};
	
	my ($bp) = (split)[2];

	my $info = (split)[7];
	
	my ($pe, $sr) = (0,0);
	if ($info =~ /PE=(\d+);/) {
		$pe = $1;
	}
	
	if ($info =~ /SR=(\d+);/) {
		$sr = $1;
	}

	my $su = $pe + $sr;
	my ($pos2) = $info =~ /END=(\d+);/;

	if ( /SVTYPE=DUP/ ){
		my ($length) = $pos2 - $pos1;
		# next if $length < 2000;
		print $delly_svs join("\t", $c1, $pos1, $pos2, $su, "type=dup") . "\n";
		$sv_count++;
	}
	
	if ( /SVTYPE=DEL/ ){
		my ($length) = $pos2 - $pos1;
		# next if $length < 2000;
		print $delly_svs join("\t", $c1, $pos1, $pos2, $su, "type=del") . "\n";
		$sv_count++;
	}
	
	if ( /SVTYPE=INV/ ){
		my ($length) = $pos2 - $pos1;
		# next if $length < 2000;
		print $delly_svs join("\t", $c1, $pos1, $pos2, $su, "type=inv") . "\n";
		$sv_count++;
	}
	
	if ( /SVTYPE=TRA/ ){
		my ($length) = $pos2 - $pos1;
		# next if $length < 2000;
		print $bnds join("\t", $c1, $pos1, $pos1, $su) . "\n";
		
		my ($c2) = $info =~ /CHR2=(.*?);/;
		
		my ($connection) = $info =~/CT=(.*?);/;
		
		
		$bp_count++;
		
		# if ($connection eq '3to3'){
# 			my $c2_stop = $pos2 + $pos1;
#
# 			if ($defined_svs and $sv_id{$bp} ){
# 				print $trans_links join("\t", $c1, $pos1, $c2, $pos2, "thickness=0.5,type=$sv_id{$bp}") . "\n";
# 			}
#
# 		}
		print $trans_links join("\t", $c1, $pos1, $pos1+1, $c2, $pos2, $pos2+1, "reads=$su") . "\n";
	
		
	}
	
	
}


say "-------------------------------------------";
say " $bp_count translocation BPs called by Delly";
say " $sv_count structural variants called by Delly";
say "-------------------------------------------";

say "Writing all Lumpy SVs to 'delly_svs.txt' [Delly SV tiles file]";
say "Writing Lumpy break points calls to 'delly_bnds.txt' [Delly BP \"heatmap\" file]";
say "Writing Lumpy translocation calls to 'delly_translocation_links.txt' [Delly translocation links]";


sub usage {
	say "Usage: delly2circos.pl [options]";
	say "  -v = delly output (VCF) [required]";
	say "  -d = user defined structural variants";
	say "  -f = filter set";
}

__DATA__
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic structural variant.">
##bcftools_concatVersion=1.3.1+htslib-1.3.1
##bcftools_concatCommand=concat -a -O v -o /data/kdi_prod/project_result/948/01.00/Analysis/Analysis/Delly/A572_73/tagged/results/A573R31.tagged.SC.gt.delly.vcf /data/kdi_prod/project_result/948/01.00/Analysis/Analysis/Delly/A572_73/tagged/results/raw_calls/A573R31.tagged.SC.DEL.gt.filt.bcf /data/kdi_prod/project_result/948/01.00/Analysis/Analysis/Delly/A572_73/tagged/results/raw_calls/A573R31.tagged.SC.DUP.gt.filt.bcf /data/kdi_prod/project_result/948/01.00/Analysis/Analysis/Delly/A572_73/tagged/results/raw_calls/A573R31.tagged.SC.INV.gt.filt.bcf /data/kdi_prod/project_result/948/01.00/Analysis/Analysis/Delly/A572_73/tagged/results/raw_calls/A573R31.tagged.SC.TRA.gt.filt.bcf
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A573R31	A573R26	A573R28	A573R30	A573R32	A573R34
2R	2944270	INV00000000	C	<INV>	0	PASS	PRECISE;SVTYPE=INV;SVMETHOD=EMBL.DELLYv0.7.6;CHR2=2R;END=16485829;PE=2;MAPQ=57;CT=3to3;CIPOS=-1,1;CIEND=-1,1;INSLEN=0;HOMLEN=0;SR=5;SRQ=0.950495;CONSENSUS=TTTAAGAAGTCCTCTTAATACATACACAAATTGCAAGAGATTAAGTAAACATGAAAGCTATCTAGTTTTATGAGATGAGATATTATGTAAAAGCTATACAC;CE=1.83874;RDRATIO=0.99711;SOMATIC	GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV	0/1:-4.17552,0,-78.0744:42:PASS:2070:7331238:2709323:5:15:2:22:3	0/0:0,-8.42652,-106.198:84:PASS:1573:7525477:2802766:5:18:0:28:0	0/0:0,-9.93158,-126.098:99:PASS:2335:9078352:3343113:5:29:0:33:0	0/0:0,-8.42723,-108.798:84:PASS:1652:7229987:2721545:5:21:0:28:0	0/0:0,-8.72785,-111.598:87:PASS:1843:7474653:2743809:5:17:0:29:0	0/0:0,-12.9418,-167.197:129:PASS:1955:9822543:3620326:5:26:0:43:0
3R	8785543	DEL00000007	A	<DEL>	.	PASS	IMPRECISE;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv0.7.6;CHR2=3R;END=8786053;PE=5;MAPQ=60;CT=3to5;CIPOS=-126,126;CIEND=-126,126;RDRATIO=0.880939;SOMATIC	GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV	0/1:-18.5609,0,-184.661:10000:PASS:197:276:73:2:33:5:0:0	0/0:0,-16.5546,-323.198:10000:PASS:192:307:96:2:55:0:0:0	0/0:0,-13.5452,-266.599:135:PASS:228:404:92:3:45:0:0:0	0/0:0,-15.3503,-302.098:154:PASS:194:318:95:2:51:0:0:0	0/0:0,-13.2435,-257.998:132:PASS:155:284:71:3:44:0:0:0	0/0:0,-20.4656,-403.996:10000:PASS:207:369:111:2:68:0:0:0
3R	31043462	DUP00000000	T	<DUP>	0	PASS	PRECISE;SVTYPE=DUP;SVMETHOD=EMBL.DELLYv0.7.6;CHR2=3R;END=31073815;PE=11;MAPQ=60;CT=5to3;CIPOS=-2,2;CIEND=-2,2;INSLEN=0;HOMLEN=1;SR=10;SRQ=1;CONSENSUS=GCCAAGGACAGTGATGAGCTCCAAGCTGAAATATTACACTGTATTGCATTTGAAAAGCAGGGTTTCGAGGAAGGTAAGTTTGGGGATCCACAGAGCGATTCTGTCCGTACTAGACTGTGTGTCAATGTC;CE=1.97828;RDRATIO=1.25974;SOMATIC	GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV	0/1:-40.8295,0,-194.527:10000:PASS:8787:23749:4816:3:14:11:56:16	0/0:0,-16.5251,-210.468:10000:PASS:9260:20081:4907:3:18:0:55:0	0/0:0,-22.271,-285.995:10000:PASS:10534:23098:5814:3:37:0:74:0	0/0:0,-13.543,-172.897:135:PASS:8672:18722:4837:3:16:0:45:0	0/0:0,-14.4455,-182.996:144:PASS:8848:18612:4623:3:14:0:48:0	0/0:0,-19.5618,-249.695:10000:PASS:11714:24835:6297:3:29:0:65:0
X	1130472	INV00000002	G	<INV>	.	PASS	IMPRECISE;SVTYPE=INV;SVMETHOD=EMBL.DELLYv0.7.6;CHR2=X;END=1131481;PE=7;MAPQ=60;CT=3to3;CIPOS=-162,162;CIEND=-162,162;RDRATIO=1.13737;SOMATIC	GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV	0/1:-30.265,0,-149.665:10000:PASS:315:502:159:2:27:7:0:0	0/0:0,-6.62264,-130.6:66:PASS:133:257:143:2:22:0:0:0	0/0:0,-6.62265,-132:66:PASS:178:291:140:2:22:0:0:0	0/0:0,-5.1175,-102:51:PASS:144:216:104:2:17:0:0:0	0/0:0,-6.01784,-116.197:60:PASS:144:287:119:2:20:0:0:0	0/0:0,-5.4182,-103.1:54:PASS:140:307:122:2:18:0:0:0
X	3136744	TRA00000000	G	<TRA>	0	PASS	PRECISE;SVTYPE=TRA;SVMETHOD=EMBL.DELLYv0.7.6;CHR2=3R;END=31063205;PE=11;MAPQ=60;CT=3to3;CIPOS=-2,2;CIEND=-2,2;INSLEN=0;HOMLEN=1;SR=10;SRQ=0.973856;CONSENSUS=CTCTCAGATTTGGAAACCCACCTTTTCGAAACAAGAAAGACCGTTTTTACCGTCGGCAAAAGGCGAGTTTTTCCAGTCCCCTGGTGAACTTGCCATGAAATCGCATTATTTGATCTAAGCCGCCGACGAGCAAAAAAGGAAATACAAAAATC;CE=1.972;RDRATIO=0;SOMATIC	GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV	0/1:-33.7829,0,-21.8831:10000:PASS:0:0:0:-1:6:10:7:10	0/0:0,-6.32031,-81.0987:63:PASS:1:1:1:1:25:0:21:0	0/0:0,-7.22315,-92.4984:72:PASS:2:2:2:1:23:0:24:0	0/0:0,-4.8154,-61.2989:48:PASS:3:3:3:1:23:0:16:0	0/0:0,-4.51443,-57.499:45:PASS:4:4:4:1:18:0:15:0	0/0:0,-6.6198,-81.8971:66:PASS:5:5:5:1:42:0:22:0
X	13954982	DEL00000010	T	<DEL>	.	PASS	IMPRECISE;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv0.7.6;CHR2=X;END=13983925;PE=4;MAPQ=45;CT=3to5;CIPOS=-242,242;CIEND=-242,242;RDRATIO=0.908944;SOMATIC	GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV	0/1:-12.9856,0,-54.1857:130:PASS:3467:6393:3551:2:10:4:0:0	0/0:0,-3.60799,-67.9956:36:PASS:3405:6380:2961:2:12:0:0:0	0/0:0,-3.61231,-68.2:36:PASS:3964:7672:3429:2:12:0:0:0	0/0:0,-4.21441,-84:42:PASS:3474:6285:2952:2:14:0:0:0	0/0:0,-2.40606,-38.5978:24:PASS:3452:6431:2766:2:8:0:0:0	0/0:0,-5.41848,-105.9:54:PASS:4508:8086:3767:2:18:0:0:0