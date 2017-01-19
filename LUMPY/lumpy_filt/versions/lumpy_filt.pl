#!/usr/bin/perl
use strict;
use warnings;

use feature qw/ say /;
use Data::Dumper;

use Getopt::Long qw/ GetOptions /;

use File::Basename;

my $vcf_file = 'test_files/tumour.normal.test.gt.vcf'; 
my $output_dir = 'test_files/';
my $help;

# Should add score threshold option
GetOptions( 'vcf=s'	        	=>		\$vcf_file, 
            'output_dir=s'      =>      \$output_dir,
            'help'              =>      \$help
	  ) or die usage();

if ($help)  { exit usage() } 

open my $in, '<', $vcf_file or die $!;
my ($name, $extention) = split(/\.([^.]+)$/, basename($vcf_file), 2);

$output_dir =~ s!/*$!/!; # Add a trailing slash

open my $hard_out, '>', $output_dir . $name . ".lumpy.hard-filter.vcf" or die $!;
open my $soft_out, '>', $output_dir . $name . ".lumpy.candidates.vcf" or die $!;

my %breaks;
while(<$in>){
    chomp;

    if ($_ =~ /^#/){   
        print $soft_out "$_\n";
        print $hard_out "$_\n";
        next;
    }
	
    my ($su) = $_ =~ /SU=(\d+);/;
    my ($id, $quality_score, $neo, $head) = (split)[2,5,9,10];
    my @neo_parts = split(/:/, $neo);
    my @head_parts = split(/:/, $head);	
	
    my ($SV_type) = $_ =~ /SVTYPE=(.*);STRANDS/;

    my ($neo_genotype, $neo_reads, $neo_pe, $neo_sr, $neo_qual, $neo_prob, $neo_li, $neo_depth) = @neo_parts[0..7];
    my ($head_genotype, $head_reads, $head_pe, $head_sr, $head_qual, $head_prob, $head_li, $head_depth) = @head_parts[0..7];

    # Rather than just skipping, it might be better to mark as "filtered" (+ colour?)
       
    next unless $head_genotype eq '0/0'; # filter out if control is NOT homozygous reference
    next unless $neo_reads >= 3;
    next unless $head_pe <= 1 and $head_sr <= 1;

    print $soft_out "$_\n";

    next unless $quality_score >= 10;
    next unless $neo_depth >= 5;
    next unless $head_pe == 0 and $head_sr == 0;

    print $hard_out "$_\n";
   
}

sub usage {
	say "********** lumpy_filter ***********";
    say "Usage: $0 [options]";
	say "  --vcf = lumpy output (see README.txt for instructions on generating this)";
	say "  --output_dir = specify name of output directory to write to";
	say "  --help";
	say "Nick Riddiford 2017";
}
