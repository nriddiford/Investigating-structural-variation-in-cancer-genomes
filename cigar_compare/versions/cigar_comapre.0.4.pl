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
my $soft_clip;

GetOptions( 'bam=s'				=>			\$bam_file,
			'debug'				=>			\$debug,
			'softclip'			=>			\$soft_clip,
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

my $debug_print;
if ($debug){
	open $debug_print, '>', $id . ".debug.log" or die $!; 
}

open my $not_filtered_bed, '>', $id . ".mapped_not_filtered.bed" or die $!;

my $count = 0;

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
		print $debug_print "_______________________________\n";
		print $debug_print "$read\n";
		print $debug_print "Droso:    $droso_cigar " . join(', ', @droso_matched) . " matched\n";
		print $debug_print (Dumper \@droso_matched);
		print $debug_print "Droso match total = $droso_match_score\n";
	}
	
	
	if (/CIGAR:/){
		my $alt_cigar_group;
		my ($alt_cigar, $alt_matched);

		# for bed output
		my ($chrom, $start, $stop, $seq) = (split)[2,3,7,9];

		my $read_length = length($seq);

		my $bed_stop = $start + $read_length;

		my $alt_mapping_flag = 0;
				
		for (split){
			if (/^[CE]D/){
				($alt_cigar_group) = $_;
				$alt_mapping_flag = 1;
			}
		}
		
		next unless $alt_mapping_flag;
		
		($alt_cigar, $alt_matched) = $alt_cigar_group =~ /CIGAR:(.*?(\d+)M.*?),/g;
				
		my @alt_matched;		
		push @alt_matched, $alt_cigar =~ /.*?(\d+)M/g;

		$alt_match_score += $_ foreach @alt_matched;
		
		if ($droso_match_score > $alt_match_score){
			print $not_filtered_bed "$chrom\t$start\t$bed_stop\n" if $debug;
			say unless $debug;
		}
		else {
			print $contamination_bed "$chrom\t$start\t$bed_stop\n";
			next;
		}
		
		if ($debug){
			print $debug_print "Alt:    $alt_cigar " . join(', ', @alt_matched) . " matched\n";
			print $debug_print (Dumper \@alt_matched);
			print $debug_print "Alt match total = $alt_match_score\n";
			($droso_match_score > $alt_match_score)
				? print $debug_print "-> Keep: $droso_match_score > $alt_match_score\n"
				: print $debug_print "-> Throw: $droso_match_score < $alt_match_score\n";
		}
		
	}
	
	else { 
		say unless $debug;
		print $debug_print "-> Keep: $droso_match_score [no alt mapping]\n" if $debug;
	}
		
	$count++;
	
	if ($debug){
		last if $count > 10000;
	}
}

sub usage {
	say "******** cigar_compare.pl ********";
    say "Usage: $0 [options]";
	say "  --bam = specify bam file processed with 'tag_reads.py'";
	say "  --debug = run in debug mode";
	say "  --test  = perform run on test dataset 'files/test.bam'";
	say "  --help = print this help message and exit\n";
	say "Nick Riddiford 2017";
}