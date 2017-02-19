#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use feature qw/ say /;

use Getopt::Long qw/ GetOptions /;

use File::Basename;

my $debug;
my $help;
my $test;
my $bam_file;

GetOptions( 'bam=s'				=>			\$bam_file,
			'debug'				=>			\$debug,
			'test'				=>			\$test, 
			'help'				=>			\$help
	  	  )	or die usage();
		  
if ($help)  { exit usage() } 

my $in;

if ($test) {
	$bam_file = 'files/test.bam';
	open $in , '<', $bam_file or die $!;
}
else {
	open $in, "samtools view -h $bam_file |" or die $!;
}

my ($name, $extention) = split(/\.([^.]+)$/, basename($bam_file), 2);

my ($id) = (split/\./, $name)[0];

open my $contamination_bed, '>', $id . ".contamination.bed" or die $!;

my $alt_count = 0;

while(<$in>){
	chomp;

	if (/^@/){
		say unless $debug;
		next;
	}

	my ($read) = (split)[0];
	my ($droso_cigar) = (split)[5];
	next unless $droso_cigar =~ /M/;
	my @droso_matched;

	push @droso_matched, $droso_cigar =~ /.*?(\d+)M/g;

	my $droso_match_score;
	my $alt_match_score;

	$droso_match_score += $_ foreach @droso_matched;

	if ($debug){
		say "_______________________________";
		say $read;
		print "Droso:    $droso_cigar " . join(', ', @droso_matched) . " matched\n";
		print Dumper \@droso_matched;
		say "Droso match total = $droso_match_score";
	}
	
	if (/CIGAR:/){
		my $alt_cigar_group;
		my ($alt_cigar, $alt_matched);
		my @split = split;

		# for bed output
		my ($chrom, $start, $stop, $seq) = (split)[2,3,7,9];

		my $read_length = length($seq);

		my $bed_stop = $start + $read_length;

		my $alt_mapping_flag = 0;
		for (@split){
			
			if (/^[CE]D/){
				$alt_cigar_group = $_;
				($alt_cigar, $alt_matched) = $alt_cigar_group =~ /CIGAR:(.*?(\d+)M.*?),/g;
				
				my @alt_matched;
			
				push @alt_matched, $alt_cigar =~ /.*?(\d+)M/g;

				$alt_match_score += $_ foreach @alt_matched;

				if ($debug){
					print "Alt:    $alt_cigar " . join(', ', @alt_matched) . " matched\n";
					print Dumper \@alt_matched;
					say "Alt match total = $alt_match_score";
				}
		
				$alt_mapping_flag = 1;
				$alt_count++;
			}
		}

		if ($alt_mapping_flag){

			if ($debug){

	 			print "Alt:    $alt_cigar    $alt_matched matched\n";
				($droso_match_score > $alt_matched) ? say "-> Keep: $droso_match_score > $alt_matched" : say "-> Throw: $droso_match_score < $alt_matched";
			}

			elsif ($droso_match_score > $alt_matched){
				say;
				print $contamination_bed "$chrom\t$start\t$bed_stop\n" unless $debug;
			}

			else { next }
		}

	 }
	 else { say unless $debug }
}

say $alt_count;

sub usage {
	say "******** cigar_compare.pl ********";
    say "Usage: $0 [options]";
	say "  --bam = specify bam file processed with 'tag_reads.py'";
	say "  --debug = run in debug mode";
	say "  --test  = run in test mode";
	say "  --help = print this help message and exit\n";
	say "Nick Riddiford 2017";
}