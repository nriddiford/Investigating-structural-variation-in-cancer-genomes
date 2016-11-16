#!/usr/bin/perl
use warnings;
use strict;

use Data::Dumper;
$Data::Dumper::Sortkeys = 1;
use feature qw(say);







while(<DATA>){
	chomp;
	next if /^#/;
	my ($info) = (split)[7];
	my @parts = split(/;/, $info);
	my ($indel, $idv, $imf, $dp, $vdp, $rpb, $mqb, $bqb, $mqsb, $sgb, $mqof) = @parts;
	
	print "***\n";
	
	print "INDEL = $indel\n" if $indel;
	print "Maximum number of reads supporting an indel = $idv\n" if $idv;
	print "Maximum fraction of reads supporting an indel = $imf\n";
	print "Raw read depth = $dp\n";
	print "Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better) = $vdp\n";
	print "Mann-Whitney U test of Read Position Bias (bigger is better) = $rpb\n";
	print "Mann-Whitney U test of Mapping Quality Bias (bigger is better) = $mqb\n";
	print "Mann-Whitney U test of Base Quality Bias (bigger is better) = $bqb\n" if $bqb;
	print "Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better) = $mqsb\n" if $mqsb;
	print "Segregation based metric = $sgb\n" if $sgb;
	print "Fraction of MQ0 reads (smaller is better) = $mqof\n" if $mqof;
	
	print "***\n";
}










__DATA__
##contig=<ID=mitochondrion_genome,length=19524>
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">
##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
##INFO=<ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">
##INFO=<ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">
##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
##INFO=<ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)">
##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=ICB,Number=1,Type=Float,Description="Inbreeding Coefficient Binomial test (bigger is better)">
##INFO=<ID=HOB,Number=1,Type=Float,Description="Bias in the number of HOMs number (smaller is better)">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">
##bcftools_callVersion=1.3.1+htslib-1.3.1
##bcftools_callCommand=call -vmO z -o /data/kdi_prod/project_result/948/01.00/Analysis/Trimmo_out/Bwa_test/Bwa_out/Mpile/HUM-1_variants.vcf.gz
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  /data/kdi_prod/project_result/948/01.00/Analysis/Trimmo_out/Bwa_test/Bwa_out/HUM-1_PE_test.bam
2L      4905    .       CAGAGAGAGAGAG   CAGAGAGAGAGAGAG 12.4322 .       INDEL;IDV=1;IMF=1;DP=1;SGB=-0.379885;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=0,0,1,0;MQ=40   GT:PL   0/1:40,3,0
2L      5974    .       C       T       9.6729  .       DP=1;SGB=-0.379885;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=0,0,0,1;MQ=60     GT:PL   0/1:37,3,0
2L      6353    .       C       T       46.3349 .       DP=2;VDB=0.18;SGB=-0.453602;MQSB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,1,1;MQ=60   GT:PL   1/1:74,6,0
2L      6631    .       A       G       9.6729  .       DP=1;SGB=-0.379885;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=0,0,0,1;MQ=60     GT:PL   0/1:37,3,0
2L      7039    .       A       T       31.39   .       DP=6;VDB=0.8;SGB=-0.453602;RPB=0.75;MQB=1;BQB=1;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=2,0,2,0;MQ=60        GT:PL   0/1:64,0,76
2L      7088    .       A       T       54      .       DP=4;VDB=0.44;SGB=-0.453602;RPB=1;MQB=1;MQSB=1;BQB=1;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=1,0,1,1;MQ=60   GT:PL   0/1:88,0,28
2L      9223    .       ATGTTGT ATGT    92      .       INDEL;IDV=3;IMF=0.6;DP=5;VDB=0.701568;SGB=-0.511536;MQSB=1;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=0,2,1,2;MQ=60     GT:PL   0/1:125,0,96
2L      10610   .       G       T       32.6366 .       DP=2;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,1,0;MQ=60   GT:PL   1/1:60,3,0
2L      13638   .       G       A       39.3353 .       DP=2;VDB=0.02;SGB=-0.453602;MQ0F=0;AC=2;AN=2;DP4=0,0,2,0;MQ=60  GT:PL   1/1:67,6,0
2L      16467   .       G       A       9.6729  .       DP=1;SGB=-0.379885;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=0,0,0,1;MQ=43     GT:PL   0/1:37,3,0
2L      16688   .       C       T       9.6729  .       DP=1;SGB=-0.379885;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=0,0,1,0;MQ=60     GT:PL   0/1:37,3,0
2L      18734   .       TATAATAATAAT    TATAATAAT       21.4353 .       INDEL;IDV=1;IMF=0.5;DP=2;VDB=0.38;SGB=-0.379885;MQSB=1;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=1,0,0,1;MQ=60 GT:PL   0/1:54,0,54
2L      19802   .       C       T       77      .       DP=3;VDB=0.0785113;SGB=-0.511536;MQSB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,1,2;MQ=60      GT:PL   1/1:105,9,0
2L      20008   .       A       G       32.6366 .       DP=2;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,0,1;MQ=60   GT:PL   1/1:60,3,0
2L      20144   .       T       A       26.1037 .       DP=4;VDB=0.32;SGB=-0.453602;RPB=1;MQB=1;MQSB=1;BQB=1;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=0,1,2,0;MQ=47   GT:PL   0/1:59,0,34
2L      20151   .       G       A       26.1037 .       DP=4;VDB=0.32;SGB=-0.453602;RPB=1;MQB=1;MQSB=1;BQB=1;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=0,1,2,0;MQ=47   GT:PL   0/1:59,0,34
2L      20242   .       A       G       9.6729  .       DP=1;SGB=-0.379885;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=0,0,0,1;MQ=52     GT:PL   0/1:37,3,0
2L      20546   .       T       C       9.6729  .       DP=1;SGB=-0.379885;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=0,0,1,0;MQ=58     GT:PL   0/1:37,3,0
2L      34143   .       C       A       9.6729  .       DP=1;SGB=-0.379885;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=0,0,0,1;MQ=60     GT:PL   0/1:37,3,0
2L      34190   .       T       G       9.6729  .       DP=1;SGB=-0.379885;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=0,0,0,1;MQ=60     GT:PL   0/1:37,3,0
2L      34345   .       A       C       21.4353 .       DP=4;SGB=-0.379885;RPB=1;MQB=1;BQB=1;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=1,0,1,0;MQ=60   GT:PL   0/1:54,0,54
2L      35460   .       C       T       9.6729  .       DP=1;SGB=-0.379885;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=0,0,1,0;MQ=60     GT:PL   0/1:37,3,0
2L      36330   .       A       T       32.6366 .       DP=2;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,0,1;MQ=60   GT:PL   1/1:60,3,0
