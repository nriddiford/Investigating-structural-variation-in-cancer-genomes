#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use feature qw/ say /;

use Bio::SeqIO;


my $vcf_file = $ARGV[0];

exit usage() unless $#ARGV == 0;

my $seqio = Bio::SeqIO->new(-file => '/Users/Nick/Documents/Curie/Data/Genomes/Dmel_v6.12/dmel-all-chromosome-r6.12.fasta', '-format' => 'Fasta');

my %chroms = qw / 2L 23513712 2R 25286936 3L 28110227 3R 32079331 4 1348131 X 23542271 Y 3667352 /;

my $genome_length;
for (keys %chroms){
	($genome_length) += $chroms{$_};
}

say $genome_length;

my %genome;
while(my $seq = $seqio->next_seq) {
	my $nucs = $seq->seq;
	my $chr = $seq->id;
	next unless $chroms{$chr}; 
	$genome{$chr} = $nucs;
}



say "Reading in VCF file: $vcf_file";
open my $VCF_in, '<', $vcf_file or die $!;

my ($name, $extention) = split(/\.([^.]+)$/, $vcf_file, 2);

open my $out, '>', $name . ".GWtrinucs.txt" or die $!;

my %snps;
my %snp_count;
my $all_snps_count = 0;
my %snp_freq;

my %genome_wide_snps;
my $sample;
say "Parsing VCF file...";
while(<$VCF_in>){
	chomp;
	next if /^##/;
	my ($chr, $pos, $ref, $alt) = (split)[0,1,3,4];
   	next if $ref eq 'N';
	my ($type) = /TYPE=(.+?)\s+/;
	$sample = (split)[9] if /^#CHROM/;

 	if ( length $ref == 1 and length $alt == 1 and $chroms{$chr} ) {
		$snps{$chr}{$pos} = [$ref, $alt];
		
		my ($trinuc) = substr($genome{$chr}, $pos - 2, 3);
		
		# say "$chr\t$nuc\t$ref=>$alt\t$pos";
		$genome_wide_snps{$trinuc}{"$ref=>$alt"}++;

		$snp_freq{$chr}{$ref}{$alt}++; # count total # transformations by type (e.g. `A` -> `G`) per chromosome
		$snp_count{$chr}++; # count total # transformations per chromosome
		$all_snps_count++; # count total # transformations
	 }
	# else { print "EXCLUDED: $chr, $ref, $alt, $type\n" }

}
say "...done";



for my $tr (sort keys %genome_wide_snps){
	for my $ref_alt (sort keys %{$genome_wide_snps{$tr}} ) {
		my ($count) = $genome_wide_snps{$tr}{$ref_alt};
		print $out "$tr\t$ref_alt\t$count\t$sample\n";
	}
}

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
	say "Usage: perl $0 <snps.vcf>";
}