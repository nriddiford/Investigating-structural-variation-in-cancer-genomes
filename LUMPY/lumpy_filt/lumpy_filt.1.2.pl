#!/usr/bin/perl
use strict;
use warnings;

use feature qw/ say /;
use Data::Dumper;

use Getopt::Long qw/ GetOptions /;

use File::Basename;

use Cwd;

my $vcf_file; 
my $output_dir = cwd();
my $help;
my $test;

# Should add score threshold option
GetOptions( 'vcf=s'	        	=>		\$vcf_file, 
            'output_dir=s'      =>      \$output_dir,
            'test'				=>		\$test,
			'help'              =>      \$help
			
	  ) or die usage();

if ($help)  { exit usage() } 

if ($test) {
	say "Running test mode...";
	$vcf_file = 'test_files/tumour.normal.test.gt.vcf';
	$output_dir = 'test_files/';
}

open my $in, '<', $vcf_file or die $!;

my ($name, $extention) = split(/\.([^.]+)$/, basename($vcf_file), 2);

$output_dir =~ s!/*$!/!; # Add a trailing slash

open my $filtered, '>', $output_dir . $name . ".lumpy.filtered.vcf" or die $!;

my $sv_count = 0;
my $pass = 0;
my $candidate_count = 0;
while(<$in>){
    chomp;

    if ($_ =~ /^#/){   
		print $filtered "$_\n";
        next;
    }
	
	my @fields = split;
    my ($su) = $_ =~ /SU=(\d+);/;
    my ($id, $quality_score, $tumour, $normal) = @fields[2,5,9,10];

    my @tumour_parts = split(/:/, $tumour);
    my @normal_parts = split(/:/, $normal);	
	
    my ($SV_type) = $_ =~ /SVTYPE=(.*);STRANDS/;

    my ($tumour_genotype, $tumour_reads, $tumour_pe, $tumour_sr, $tumour_qual, $tumour_prob, $tumour_li, $tumour_depth) = @tumour_parts[0..7];
    my ($normal_genotype, $normal_reads, $normal_pe, $normal_sr, $normal_qual, $normal_prob, $normal_li, $normal_depth) = @normal_parts[0..7];
	       
    next unless $normal_genotype eq '0/0'; 			# ignore if normal sample is NOT homozygous reference
    next unless $tumour_reads >= 3;					# ignore if fewer than 3 reads in tumour sample supporting SV
    next unless $normal_pe < 2;						# ignore if more than 1 discordant read in normal sample
	next unless $normal_sr < 2;						# ignore if more than 1 split read in normal sample
	
	my @filter_reasons;
		
	if ($quality_score < 10 ){
    	push @filter_reasons, 'qual=' . $quality_score;
		$pass = 1;
    }
	
  	if ($tumour_depth < 5 ){
		push @filter_reasons, 'depth=' . $tumour_depth;
		$pass = 1;
	}	
	
	if ($normal_pe > 0 ){
		push @filter_reasons, 'CPE=' . $normal_pe;
		$pass = 1;
	}
	
	if ($normal_sr > 0 ) {
		push @filter_reasons, 'CSR=' . $normal_sr;
		$pass = 1;
	}
			
	elsif ( 
		 $quality_score	>=	10	and
		 $tumour_depth	>=	5	and
		 $normal_reads	==	0
	   ) {
			
	 	print $filtered join("\t", @fields[0..5], "PASS", @fields[7..10]) . "\n";
		$sv_count++;
	}
	
	$candidate_count++ if $pass;
	
	next unless scalar @filter_reasons > 0;
	
	my @filter_print = join(",", @filter_reasons);

	print $filtered join("\t", @fields[0..5], @filter_print, @fields[7..10]) . "\n";
 
}

say "$sv_count structual variants passed hard filter";
my $candidates = $candidate_count - $sv_count;
say "$candidates candidate structual variants";

sub usage {
	say "********** lumpy_filter ***********";
    say "Usage: $0 [options]";
	say "  --vcf = lumpy output (see README.txt for instructions on generating this)";
	say "  --output_dir = specify name of output directory to write to";
	say "  --test = run on test data";
	say "  --help\n";
	say "Nick Riddiford 2017";
}
