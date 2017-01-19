#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use feature qw/ say /;

use Bio::SeqIO;

use Getopt::Long qw/ GetOptions /;

my $genome_file = '/Users/Nick/Documents/Curie/Data/Genomes/Dmel_v6.12/Dmel_6.12.fasta';
my $vcf_file = '/Users/Nick/Desktop/Analysis/Freebayes/HUM-1.TEx.q4.frba.vcf';
my $chrom_out_file = 'chroms.trinucs.txt';
my $genome_out_file = 'GW.trinucs.txt';
my $debug;
my $quiet;
my $help;


# Should add score threshold option
GetOptions( 'genome=s'     			=>    	\$genome_file,
			'vcf=s'					=>		\$vcf_file,
		   	'all-snps=s'    		=>	 	\$genome_out_file,
		   	'chrom-snps=s'    		=>	 	\$chrom_out_file,
		   	'help'         			=>   	\$help,
		   	'quiet'        			=>   	\$quiet,
           	'debug'        			=>    	\$debug
) or die usage();

if ($help)  { exit usage() } 
if ($quiet) { say "Running in quiet mode" }
if ($debug) { say "Running in debug mode" }


open my $chrom_out, '>>',  $chrom_out_file or die $!;
open my $genome_out, '>>',  $genome_out_file or die $!;

my $seqio = Bio::SeqIO->new(-file => "$genome_file", '-format' => 'Fasta');

my %chroms = qw / 2L 23513712 2R 25286936 3L 28110227 3R 32079331 4 1348131 X 23542271 Y 3667352 /;


my $genome_length;

for (keys %chroms){
	($genome_length) += $chroms{$_};
}

say "Reading in genome: $genome_file";

my %genome;

while(my $seq = $seqio->next_seq) {
	my $nucs = $seq->seq;
	my $chr = $seq->id;
	next unless $chroms{$chr}; 
	$genome{$chr} = $nucs;
}

my ($name, $extention) = split(/\.([^.]+)$/, $vcf_file, 2);

my %snps;
my %snp_count;
my $all_snps_count = 0;
my %snp_freq;

my %genome_wide_snps;
my $sample;

my %snps_by_chrom;
my %tri_count;

say "Reading in VCF file: $vcf_file";
open my $VCF_in, '<', $vcf_file or die $!;

say "Parsing VCF file...";

while(<$VCF_in>){
	chomp;
	next if /^##/;
	my ($chr, $pos, $ref, $alt, $quality, $info) = (split)[0,1,3,4,5,7];
   	next if $ref eq 'N';
	my ($type) = /TYPE=(.+?)\s+/;
	$sample = (split)[9] if /^#CHROM/;
    
	my ($read_depth) = $info =~ /DP=(\d+);/;
	my ($mmq, $mmq_ref) = $info =~ /MQM=(\d+.?\d*);MQMR=(\d+.?\d*);/;
	
 	if ( length $ref == 1 and length $alt == 1 and $chroms{$chr} ) {
				
		# Hard filters (see: 
		# http://gatkforums.broadinstitute.org/gatk/discussion/1186/best-practice-variant-detection-with-the-gatk-v4-for-release-2-0
		# https://www.biostars.org/p/85400/
		# https://bcbio.wordpress.com/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/
		# QUAL < 30
		# MQ < 40.0
		# DP < 10
		
		next if $quality < 30;
		next if $read_depth < 10;
		next if $mmq < 40;
		
		$snps{$chr}{$pos} = [$ref, $alt];
		
		my ($trinuc) = substr( $genome{$chr}, $pos - 2, 3 );
	
		if ($trinuc =~ /N/){
			say "excluding $trinuc";
			next;
		}
				
		$genome_wide_snps{$trinuc}{"$ref>$alt"}++;		# count genome-wide trinuc transformations: "A>C"

		$snps_by_chrom{$chr}{$trinuc}{"$ref>$alt"}++;	# count chromosome-wide trinuc transformations: 2L "A>C"
		
		$snp_freq{$chr}{$ref}{$alt}++; 					# count total transformations by type (e.g. `A` -> `G`) per chromosome
		$snp_count{$chr}++; 							# count total transformations per chromosome
		$all_snps_count++; 								# count total transformations
		$tri_count{$chr}++;								# count total trinucs per chromosome
		
		my $snp_count = $snp_freq{$chr}{$ref}{$alt};
		my ($mut_cont) = eval sprintf('%.1f', $snp_count/$snp_count{$chr} * 100);
		
		if ($debug){
			say "SNP:                            $ref";
			say "Position:			$chr\:$pos";
			say "Quality score:			$quality";
			say "Mapping quality:		$mmq";
			say "Transition:                     $ref>$alt";
			say "Contribution to mutation load:  $mut_cont\%";
			say "Trinucleotide:                  $trinuc";
			say "***********";
		}
		
	 }

}
say "...done";

say "Printing out SNPS per chromosome to '$chrom_out_file'...";

for my $chr (sort keys %snps_by_chrom){
	for my $tri (sort keys %{$snps_by_chrom{$chr}} ) {
		for my $ref_alt (sort keys %{$snps_by_chrom{$chr}{$tri}} ){
			my ($count) = $snps_by_chrom{$chr}{$tri}{$ref_alt};
			my ($freq) = eval sprintf('%.2f', ( $count/$snp_count{$chr} ) * 100);
			print $chrom_out join("\t", $chr, $tri, $ref_alt, $freq, $sample) . "\n";
		}
			
	}
}
say "...done";

say "Printing out genome-wide SNPS '$genome_out_file'...";

for my $tr (sort keys %genome_wide_snps){
	for my $ref_alt (sort keys %{$genome_wide_snps{$tr}} ) {
		my ($count) = $genome_wide_snps{$tr}{$ref_alt};
		my ($freq) = eval sprintf('%.2f', ( $count/$all_snps_count ) * 100);
		print $genome_out join("\t", $tr, $ref_alt, $freq, $sample) . "\n";
	}
}
say "...done";


# say "Finding mutations in index and calculating frequencies";
# print "\nFound $all_snps_count snps across whole genome\n\n";
#
# my %tri;
# for my $chromosome ( sort keys %snps ){
#  	say "Found $snp_count{$chromosome} snps on $chromosome (which is $chroms{$chromosome} bases)";
#  	my ($freq) = $snp_count{$chromosome}/$chroms{$chromosome};
#  	say "Mutation rate on chromosome:    $freq";
#          for my $position ( sort keys %{$snps{$chromosome}} ){
# 		 	my $reference_from_index = $index{$chromosome}{$position};
# 			my $reference_from_vcf = @{$snps{$chromosome}{$position}}[0];
#
# 			my $alt_from_vcf = @{$snps{$chromosome}{$position}}[1];
#
# 			unless ($reference_from_index eq $reference_from_vcf ){
# 				say "$reference_from_index ne $reference_from_vcf!!";
# 				next;
# 			}
#
# 			my $before = $index{$chromosome}{$position - 1};
# 			my $after = $index{$chromosome}{$position + 1};
#
# 			my $trinucleotide = $before . $reference_from_index . $after;
#
# 			# store freqs in hash for later
# 			$tri{$chromosome}{$trinucleotide}++;
#
#  			my $count = $snp_freq{$chromosome}{$reference_from_vcf}{$alt_from_vcf};
# 			my $percentage = eval sprintf('%.1f', $count/$snp_count{$chromosome} * 100);
#
# 			say "SNP:                            $reference_from_vcf";
# 			say "Transition:                     $reference_from_vcf>$alt_from_vcf";
# 			say "Contribution to mutation load:  $percentage\%";
# 			say "Trinucleotide:                  $trinucleotide";
# 			say "***********";
# 		}
#  	print "\n";
# }
# say "...done";
#
# say "writing trinucleotides to $out...";
# for my $ch (sort keys %tri){
# 	for my $seq (sort keys %{$tri{$ch}}) {
# 		my $count = $tri{$ch}{$seq};
# 		my $percentage = eval sprintf('%.1f', $count/$snp_count{$ch} * 100);
# 		print $out "$ch $seq $percentage\n";
# 	}
# }
#
# say "...done";

sub usage {
    say "********** compare_snps ***********";
    say "Usage: $0 [options]";
	say "--genome = genome fasta file";
	say "--all-snps = specify name of output file for all snps";
	say "--chrom-snps = specify name of output file for snps per chromosome";
	say "--help";
	say "--quiet";
	say "--debug";
	say "Nick Riddiford 2016";
}