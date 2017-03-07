#!/usr/bin/perl

use strict;
use warnings;

use feature qw/ say /;

use Data::Dumper;

use Getopt::Long qw/ GetOptions /;

my $characterised_svs;
my $vcf;
my $filter;

GetOptions( 'vcf_in=s'					=>			\$vcf,
			'defined_variants=s'		=>			\$characterised_svs,
			'filter_set=s'				=>			\$filter
	  	  )	or die usage();
		  
exit usage() unless $vcf;

my $characterised_varaints;

my %defined;

if ($characterised_svs) {
	open $characterised_varaints , '<', $characterised_svs or die $!;
	
	while(<$characterised_varaints>){
		chomp;
		my ($id, $type, $flag) = split;
		my @flags = split(',', $flag);
		$defined{$id} = [$type, @flags];
	}
	
	say "\nAnnotated:";
	print Dumper \%defined;
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

open my $trans_links, '>', 'lumpy_translocation_links.txt' or die $!;
open my $bnds, '>', 'lumpy_bnds.txt' or die $!;

open my $lumpy_svs, '>', 'lumpy_svs.txt' or die $!;

open my $circos_svs, '>', 'svs_circos.vcf' or die $!;


my @keys = qw / 2L 2R 3L 3R 4 X Y /;
my %chrom_filt;
$chrom_filt{$_} = 1 for (@keys);


my $sv_count = 0;
my $bp_count = 0;
while(<$in>){
	chomp;
	
	if (/^#/){
		print $circos_svs "$_\n";
		next;
	}
	
	my ($chrom1, $start, $pass) = (split)[0,1,6];
	next unless $chrom_filt{$chrom1};
	
	# Get the SV ID
	my ($bp) = (split)[2];
	my ($id) = $bp =~ /(\d+)_?/;
	
	unless ($pass eq 'PASS' or $pass eq 'CPE=1' or $pass eq 'CSR=1' or $defined{$id}){
		next;
	}
	
	# If it's in the filter list print the ID and reason and skip
	if ( $id_filter{$id} ){
		next;
	}
	
	my $info = (split)[7];
	my ($type) = $info =~ /SVTYPE=(.*?);/;
	
	my @flags;
	
	
	if ( $defined{$id} ){
			$type = @{$defined{$id}}[0];
			push @flags, @{$defined{$id}}[1];
	}
	
	my ($su) = $info =~ /SU=(\d+);/;
	push @flags, "stroke_thickness=0" if $pass ne 'PASS';
	
	my $stop;
	
	if ( $type eq 'DUP' ){
		push @flags, "type=dup";
	}
	
	if ( $type eq 'DEL' ){
		push @flags, "type=del";
	}
	if ( $type eq 'INV' ){
		push @flags, "type=inv";
	}
	
	my $flags_for_circos = join(",", @flags );
	
	unless ($type eq 'BND'){
		my ($length, $stop, $chrom2);
		if ( $defined{$id} ){
			my ($alt) = (split)[4];
			$alt =~ s/[\[\]N]//g;
			($chrom2, $stop) = $alt =~ /(.+)\:(\d+)/;
			$length = $stop - $start;
		}
		else {
			($length) = $info =~ /;SVLEN=-?(\d+);/;
			($stop) = $info =~ /END=(\d+);/;
		}
		# next if $length < 2000;
		print $lumpy_svs join("\t", $chrom1, $start, $stop, $su, $flags_for_circos) . "\n";
		print $circos_svs "$_\n";
		$sv_count++;
	}
	
	
	# if ( $type eq 'DUP' ){
	# 	push @flags, "type=dup";
	# 	my $
	#
	# 	# $pass eq 'PASS' ? print $lumpy_svs join("\t", $chrom1, $start, $stop, $su, @flags) . "\n" : print $lumpy_svs join("\t", $chrom1, $start, $stop, $su, "stroke_thickness=0,type=dup") . "\n" ;
	#
	# }
	#
	# if ( $type eq 'DEL' ){
	# 	push @flags, "type=del";
	# 	# my ($length) = $info =~ /;SVLEN=-?(\d+);/;
	# 	# next if $length < 2000;
	#
	# 	$pass eq 'PASS' ? print $lumpy_svs join("\t", $chrom1, $start, $stop, $su, "type=del") . "\n" : print $lumpy_svs join("\t", $chrom1, $start, $stop, $su, "stroke_thickness=0,type=del") . "\n" ;
	# 	print $lumpy_svs join("\t", $chrom1, $start, $stop, $su, @flags) . "\n";
	# }
	#
	# if ( $type eq 'INV' ){
	# 	my ($length) = $info =~ /;SVLEN=-?(\d+);/;
	# 	next if $length < 2000;
	# 	$pass eq 'PASS' ? print $lumpy_svs join("\t", $chrom1, $start, $stop, $su, "type=inv") . "\n" : print $lumpy_svs join("\t", $chrom1, $start, $stop, $su, "stroke_thickness=0,type=inv") . "\n" ;
	# }
	

	if ( $type eq 'BND'){
		$type = 'TRA';
		unless ($pass eq 'PASS' or $defined{$id}) {
			next;
		}
		push @flags, "type=tra";
		
		my $bp1 = $start;
		# Get translocation bp
		my ($alt) = (split)[4];
		$alt =~ s/[\[\]N]//g;
		
		my ($chrom2, $bp2) = $alt =~ /(.+)\:(\d+)/;
		# skip unless transloction is to a whole chromosome
		next if not $chrom_filt{$chrom2};
		
		if ($start < $bp2){
			# next if ($bp2 - $bp1) < 500;
		}
		elsif ($start > $bp2){
			# next if ($bp1 - $bp2) < 500;
		}
		
		my $start_plus = $bp1 + 1;
		my $stop_plus = $bp2 + 1;
		
		if ($su <= 4){
			push @flags, "thickness=0.5";
			my $flags_for_trans = join(",", @flags );
			print $trans_links join("\t", $chrom1, $bp1, $start_plus, $chrom2, $bp2, $stop_plus, $flags_for_trans) . "\n";
			# if ($characterised_svs and $defined{$id} ){
# 				print $trans_links join("\t", $chrom1, $bp1, $start_plus, $chrom2, $bp2, $stop_plus, $flags_for_circos) . "\n";
# 			}
# 			else { print $trans_links join("\t", $chrom1, $bp1, $start_plus, $chrom2, $bp2, $stop_plus, "thickness=0.5") . "\n"}
 		}
		
		elsif ($su >= 5 ){
			push @flags, "thickness=4";
			my $flags_for_trans = join(",", @flags );
			print $trans_links join("\t", $chrom1, $bp1, $start_plus, $chrom2, $bp2, $stop_plus, $flags_for_trans) . "\n";
			
			# if ($characterised_svs and $defined{$id} ){
# 				print $trans_links join("\t", $chrom1, $bp1, $start_plus, $chrom2, $bp2, $stop_plus, "thickness=4,type=$defined{$bp}") . "\n";
# 			}
# 			else { print $trans_links join("\t", $chrom1, $bp1, $start_plus, $chrom2, $bp2, $stop_plus, "thickness=4") . "\n"}
		}
		my $bp1_end = $bp1 += 1;
		print $bnds join("\t", $chrom1, $bp1, $bp1, $su) . "\n";
		$bp_count++;
		print $circos_svs "$_\n";
	}
	
}

say "-------------------------------------------";
say " $bp_count translocation BPs called by Lumpy";
say " $sv_count structural variants called by Lumpy";
say "-------------------------------------------";

say "Writing all Lumpy SVs to 'lumpy_svs.txt' [Lumpy SV tiles file]";
say "Writing Lumpy break points calls to 'lumpy_bnds.txt' [Lumpy BP \"heatmap\" file]";
say "Writing Lumpy translocation calls to 'lumpy_translocation_links.txt' [Lumpy translocation links]";

sub usage {
	say "Usage: lumpy2circos.pl [options]";
	say "  -v = lumpy output (VCF) [required]";
	say "  -d = user characterised structural variants. This must be a two column file: Lumpy_ID	Type [e.g. '123	INV']";
	say "  -f = user defined exclude list. This must be a two column file: Lumpy_ID	Reason [e.g. '123	False positve']";
}

__DATA__
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HUM-1	HUM-3
2L	2704378	133	N	<DEL>	0.00	qual=0.00	SVTYPE=DEL;SVLEN=-917655;END=3622033;STRANDS=+-:5;IMPRECISE;CIPOS=-7,553;CIEND=-566,5;CIPOS95=0,73;CIEND95=-73,0;SU=5;PE=5;SR=0	GT:SU:PE:SR:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB	0/0:5:5:0:96:0.00:-0,-10,-32:34:33:0:32:0:19:0:0:13:0:0	0/0:0:0:0:200:0.00:-0,-26,-85:86:86:0:85:0:54:0:0:31:0:0
2L	803778	173_1	N	N]2L:4334462]	0.00	qual=0.00	SVTYPE=BND;STRANDS=++:4;IMPRECISE;CIPOS=-10,9;CIEND=-10,9;CIPOS95=-1,1;CIEND95=-1,1;MATEID=173_2;EVENT=173;SU=4;PE=2;SR=2	GT:SU:PE:SR:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB	0/0:4:2:2:45:0.00:-3,-8,-30:38:34:3:33:2:21:1:0:12:1:0.057	0/0:0:0:0:200:0.00:-0,-24,-81:82:82:0:81:0:49:0:0:32:0:0
2L	4334462	173_2	N	N]2L:803778]	0.00	qual=0.00	SVTYPE=BND;STRANDS=++:4;IMPRECISE;CIPOS=-10,9;CIEND=-10,9;CIPOS95=-1,1;CIEND95=-1,1;MATEID=173_1;EVENT=173;SECONDARY;SU=4;PE=2;SR=2	GT:SU:PE:SR:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB	0/0:4:2:2:45:0.00:-3,-8,-30:38:34:3:33:2:21:1:0:12:1:0.057	0/0:0:0:0:200:0.00:-0,-24,-81:82:82:0:81:0:49:0:0:32:0:0
2L	1143324	207_1	N	[2L:5461159[N	0.00	qual=0.00	SVTYPE=BND;STRANDS=--:4;IMPRECISE;CIPOS=-590,9;CIEND=-590,9;CIPOS95=-94,1;CIEND95=-94,1;MATEID=207_2;EVENT=207;SU=4;PE=4;SR=0	GT:SU:PE:SR:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB	0/0:4:4:0:54:0.00:-0,-5,-18:20:19:0:18:0:8:0:0:10:0:0	0/0:0:0:0:159:0.00:-0,-16,-54:55:55:0:54:0:22:0:0:32:0:0
2L	5461159	207_2	N	[2L:1143324[N	0.00	PASS	SVTYPE=BND;STRANDS=--:4;IMPRECISE;CIPOS=-590,9;CIEND=-590,9;CIPOS95=-94,1;CIEND95=-94,1;MATEID=207_1;EVENT=207;SECONDARY;SU=4;PE=4;SR=0	GT:SU:PE:SR:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB	0/0:4:4:0:54:0.00:-0,-5,-18:20:19:0:18:0:8:0:0:10:0:0	0/0:0:0:0:159:0.00:-0,-16,-54:55:55:0:54:0:22:0:0:32:0:0
2L	5082362	213_1	N	N]2L:5569537]	0.00	qual=0.00	SVTYPE=BND;STRANDS=++:4;IMPRECISE;CIPOS=-10,573;CIEND=-9,583;CIPOS95=-1,94;CIEND95=0,94;MATEID=213_2;EVENT=213;SU=4;PE=4;SR=0	GT:SU:PE:SR:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB	0/0:4:4:0:114:0.00:-0,-11,-38:40:39:0:38:0:25:0:0:13:0:0	0/0:0:0:0:200:0.00:-0,-24,-79:80:80:0:79:0:48:0:0:31:0:0
2L	5569537	213_2	N	N]2L:5082362]	0.00	qual=0.00	SVTYPE=BND;STRANDS=++:4;IMPRECISE;CIPOS=-9,583;CIEND=-10,573;CIPOS95=0,94;CIEND95=-1,94;MATEID=213_1;EVENT=213;SECONDARY;SU=4;PE=4;SR=0	GT:SU:PE:SR:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB	0/0:4:4:0:114:0.00:-0,-11,-38:40:39:0:38:0:25:0:0:13:0:0	0/0:0:0:0:200:0.00:-0,-24,-79:80:80:0:79:0:48:0:0:31:0:0
2L	1121578	235_1	N	[2L:6085249[N	0.00	qual=0.00	SVTYPE=BND;STRANDS=--:4;IMPRECISE;CIPOS=-4,9;CIEND=-10,7;CIPOS95=0,1;CIEND95=-2,0;MATEID=235_2;EVENT=235;SU=4;PE=0;SR=4	GT:SU:PE:SR:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB	0/0:4:0:4:120:0.00:-0,-12,-40:42:41:1:40:0:25:0:0:15:0:0	0/0:0:0:0:200:0.00:-0,-27,-91:92:92:0:91:0:55:0:0:36:0:0
2L	6085249	235_2	N	[2L:1121578[N	0.00	PASS	SVTYPE=BND;STRANDS=--:4;IMPRECISE;CIPOS=-10,7;CIEND=-4,9;CIPOS95=-2,0;CIEND95=0,1;MATEID=235_1;EVENT=235;SECONDARY;SU=4;PE=0;SR=4	GT:SU:PE:SR:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB	0/0:4:0:4:120:0.00:-0,-12,-40:42:41:1:40:0:25:0:0:15:0:0	0/0:0:0:0:200:0.00:-0,-27,-91:92:92:0:91:0:55:0:0:36:0:0
2L	7928988	348_1	N	[2L:8630730[N	0.00	qual=0.00	SVTYPE=BND;STRANDS=--:4;IMPRECISE;CIPOS=-2,9;CIEND=-9,8;CIPOS95=0,8;CIEND95=0,0;MATEID=348_2;EVENT=348;SU=4;PE=0;SR=4	GT:SU:PE:SR:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB	0/0:4:0:4:117:0.00:-0,-12,-39:42:40:1:39:0:26:0:0:13:0:0	0/0:0:0:0:200:0.00:-0,-21,-70:71:71:0:70:0:36:0:0:34:0:0
2L	10280205	614	N	<DUP>	1.15	qual=1.15	SVTYPE=DUP;SVLEN=4448130;END=14728335;STRANDS=-+:4;IMPRECISE;CIPOS=-10,9;CIEND=-10,9;CIPOS95=-1,0;CIEND95=-1,1;SU=4;PE=2;SR=2	GT:SU:PE:SR:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB	0/0:4:2:2:6:1.15:-1,-1,-3:14:12:1:12:1:0:0:0:12:1:0.077	0/0:0:0:0:140:0.00:-0,-14,-25:83:83:0:82:0:49:0:0:33:0:0