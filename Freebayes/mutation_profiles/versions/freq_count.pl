#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use feature qw/ say /;

my $vcf_file = $ARGV[0];
my $index_file = $ARGV[1];

exit usage() unless $#ARGV == 1;

say "Running new version";

say "Reading in VCF file: $vcf_file";
open my $VCF_in, '<', $vcf_file or die $!;

my ($name, $extention) = split(/\.([^.]+)$/, $vcf_file, 2);

say "Reading in genome index: $index_file";
open my $index, '<', $index_file or die $!;

# my %chroms = qw / X 23542271 /;
my %chroms = qw / 2L 23513712 2R 25286936 3L 28110227 3R 32079331 4 1348131 X 23542271 Y 3667352 /;

my ($number_of_lines) = `wc -l $index_file | tr -s ' ' | cut -d ' ' -f 2`;
chomp($number_of_lines);

my @completion = qw/ 10 20 30 40 50 60 70 80 90 100 /;

say "Building hash table for genome...";
my %index;
while(<$index>){
	chomp;
	my ($chr, $nuc, $pos) = split;
	if ($chroms{$chr}){
		$index{$chr}{$pos} = $nuc;
	}
	my ($line_no) = $.;
	my ($completed) = ($line_no / $number_of_lines) * 100;
	for my $perc (@completion){
		say "$perc\%..." if ($completed == $perc);
	}
}


say "...done";

open my $out, '>', $name . ".data.txt" or die $!;

my %snps;
my %snp_count;
my $all_snps_count = 0;
my %snp_freq;
# my %sort_by_count;


say "Parsing VCF file...";
while(<$VCF_in>){
	chomp;
	next if /^#/;
	# my ($chr, $ref, $alt) = (split)[0,3,4];
	my ($chr, $pos, $ref, $alt) = (split)[0,1,3,4];
   	next if $ref eq 'N';
	my ($type) = /TYPE=(.+?)\s+/;

 	if ( length $ref == 1 and length $alt == 1 and $chroms{$chr} ) {
		$snps{$chr}{$pos} = [$ref, $alt];

		# $snps{$chr}{$ref}{$alt}++;
		# $sort_by_count{$chr} = $snps{$chr}{$ref}{$alt};

		$snp_freq{$chr}{$ref}{$alt}++; # count total # transformations by type (e.g. `A` -> `G`) per chromosome

		$snp_count{$chr}++; # count total # transformations per chromosome
 		$all_snps_count++; # count total # transformations
	}
	# else { print "EXCLUDED: $chr, $ref, $alt, $type\n" }

}
say "...done";

say "Finding mutations in index and calculating frequencies";
print "\nFound $all_snps_count snps across whole genome\n\n";

my %tri;
for my $chromosome ( sort keys %snps ){
 	say "Found $snp_count{$chromosome} snps on $chromosome (which is $chroms{$chromosome} bases)";
 	my ($freq) = $snp_count{$chromosome}/$chroms{$chromosome};
 	say "Mutation rate on chromosome:    $freq";
         for my $position ( sort keys %{$snps{$chromosome}} ){
		 	my $reference_from_index = $index{$chromosome}{$position};
			my $reference_from_vcf = @{$snps{$chromosome}{$position}}[0];

			my $alt_from_vcf = @{$snps{$chromosome}{$position}}[1];

			unless ($reference_from_index eq $reference_from_vcf ){
				say "$reference_from_index ne $reference_from_vcf!!";
				next;
			}

			my $before = $index{$chromosome}{$position - 1};
			my $after = $index{$chromosome}{$position + 1};

			my $trinucleotide = $before . $reference_from_index . $after;

			# store freqs in hash for later
			$tri{$chromosome}{$trinucleotide}++;

 			my $count = $snp_freq{$chromosome}{$reference_from_vcf}{$alt_from_vcf};
			my $percentage = eval sprintf('%.1f', $count/$snp_count{$chromosome} * 100);

			say "SNP:                            $reference_from_vcf";
			say "Transition:                     $reference_from_vcf>$alt_from_vcf";
			say "Contribution to mutation load:  $percentage\%";
			say "Trinucleotide:                  $trinucleotide";
			say "***********";
		}
 	print "\n";
}
say "...done";

say "writing trinucleotides to $out...";
for my $ch (sort keys %tri){
	for my $seq (sort keys %{$tri{$ch}}) {
		my $count = $tri{$ch}{$seq};
		my $percentage = eval sprintf('%.1f', $count/$snp_count{$ch} * 100);
		print $out "$ch $seq $percentage\n";
	}
}

say "...done";

sub usage {
	say "Usage: perl $0 <snps.vcf> <genome.index";
}