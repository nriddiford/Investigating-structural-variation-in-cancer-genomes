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
			'softclip'			=>			\$soft_clip,
			'test'				=>			\$test, 
			'debug'				=>			\$debug,
			'help'				=>			\$help
	  	  )	or die usage();
		  
if ($help)  { exit usage() } 
exit usage() unless $bam_file or $test;

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
open my $te_bed, '>', $id . ".TE.bed" or die $!;

my ($debug_print, $sc_filter_print, $not_filtered_bed);
if ($debug){
	open $debug_print, '>', $id . ".debug.txt" or die $!; 
	open $sc_filter_print, '>', $id . ".softclip_filtered.bed" or die $!;
	open $not_filtered_bed, '>', $id . ".mapped_not_filtered.bed" or die $!;
}


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
	
	$count++;

	# For bed output
	my ($chrom, $start, $stop, $seq) = (split)[2,3,7,9];
	my $read_length = length($seq);
	my $bed_stop = $start + $read_length;
	
	my ($droso_match_score, @droso_matched) = sum_matches($droso_cigar);
	
	debug($read, $droso_cigar, $droso_match_score, "Droso", @droso_matched) if $debug;
	
	# First, skip any reads exceeding the soft clippoing threshold of SC > 15 M < 30 SC > 15 (if run in soft clip mode '-s')
	if ($soft_clip and $droso_cigar =~ /\d+S\d+M\d+S/){	
		my ($sc_1, $m, $sc_2) = $droso_cigar =~ /(\d+)S(\d+)M(\d+)S/;
			
			if ($sc_1 > 15 and $m < 30 and $sc_2 > 15){
				# Make a bed file for sc reads removed if running in debug mode
				print $sc_filter_print "$chrom\t$start\t$bed_stop\n" if $debug;
				
				# Say why this read is being skipped
				print $debug_print "-> Throw read: $droso_cigar [exceeded double-ended soft-clipping threshold]\n" if $debug;
				next;
			}
	}
	
	# Next, if the read line contains "CIGAR" it's been modified by 'read_tagger.py' and tagged as mapping to an alt geneome 
	
	if (/CIGAR:/){
		my $alt_cigar_group;
		my $non_filter_alt_cigar_group;
		my $mate_cigar_group;

		my $alt_mapping_flag = 0;
		my $non_filter_set = 0;
		my $mate_pass = 0;
		
		# Break the read line up and check to see if each block starts with a tag descriptor (e.g. 'CD')
		
			for (split){
			
				# Primary tag group (Read = A, mate = B). In this instance, this is the 'non-filter' set
				if (/^[A]D/){
					$non_filter_alt_cigar_group = $_;
					$non_filter_set = 1;
				}
			
				# Check secondary tag groups 
				elsif (/^[CEG]D/){
					($alt_cigar_group) = $_;
					 $alt_mapping_flag = 1;
				 }
			
				 # We want to keep the mate (who's also been tagged)
				elsif (/^[BDFH]D/){
					($mate_cigar_group) = $_;
					$mate_pass = 1;
				}
			
			}
						
			# For the non-filter set, we want to keep the tagged reads, and make a bed file showing their locations
			if ($non_filter_set){

				my ($non_filt_alt_cigar) = $non_filter_alt_cigar_group =~ /CIGAR:(.*?\d+M.*?),/g;
	
				my ($non_filt_alt_matched_score, @non_filt_matched) = sum_matches($non_filt_alt_cigar);
								
				debug($read, $non_filt_alt_cigar, $non_filt_alt_matched_score, "TE", @non_filt_matched) if $debug;
				
				say unless $debug;
				
				print $te_bed "$chrom\t$start\t$bed_stop\n";
			}
		
			# next unless $alt_mapping_flag;
			
			# For the alt mappers, we want to check that they map better to the alt reference than to drosophila and decide on whether to discard or not
			if ($alt_mapping_flag){
				
				my ($alt_cigar) = $alt_cigar_group =~ /CIGAR:(.*?\d+M.*?),/g;
				
				my ($alt_match_score, @alt_matched) = sum_matches($alt_cigar);
															
				# If the total map score is higer in droso than in the alt ref
				if ($droso_match_score > $alt_match_score){
					print $not_filtered_bed "$chrom\t$start\t$bed_stop\n" if $debug;
					debug($read, $alt_cigar, $alt_match_score, "false_alt", @alt_matched) if $debug;
					say unless $debug;
				}
				
				# If not, skip read and make bed showing contamination
				elsif ($droso_match_score < $alt_match_score){
					print $contamination_bed "$chrom\t$start\t$bed_stop\n";
					debug($read, $alt_cigar, $alt_match_score, "true_alt", @alt_matched) if $debug;
					next;
				}
				
			}
		
			elsif ($mate_pass){
				my ($mate_cigar) = $mate_cigar_group =~ /CIGAR:(.*?\d+M.*?),/g;
				my ($mate_match_score, @mate_match) = sum_matches($mate_cigar);
				say unless $debug;
				
				debug($read, $mate_cigar, $mate_match_score, "mate", @mate_match) if $debug;
				next;
			}
		}
	
	# Pass all reads that aren't tagged
	else {
		print $debug_print "-> Keep read: [no alt mapping]\n" if $debug;
		say unless $debug;
	}
			
	if ($debug){
		last if $count == 10000;
	}
	
}

sub debug {
	my ($read, $cigar, $score, $type, @matched) = @_;
	
	if ($type eq 'Droso'){
		print $debug_print "-----------------------------------------------\n";
		print $debug_print "$read\n";
	}
	
	print $debug_print "$type cigar: $cigar\n" . "Matches to $type: " . join(', ', @matched) . "\n";
	print $debug_print "  " . (Dumper \@matched);
	print $debug_print "$type match total = $score\n";
	print $debug_print "------------\n";
	
	if ($type eq 'TE'){
		print $debug_print "This read has been tagged as mapping to a TE\n";
		print $debug_print "-> Keep read\n";
	}
	if ($type eq 'false_alt'){
		print $debug_print "This read has been tagged as mapping to a bacterial geneome\n";
		print $debug_print "-> Keep read: [Maps better to Drosophila]\n"
	}
	if ($type eq 'true_alt'){
		print $debug_print "This read has been tagged as mapping to a bacterial geneome\n";
		print $debug_print "-> Throw read: [Contamination]\n";
	} 
	if ($type eq 'mate'){
		print $debug_print "This read has been tagged becuase its mate maps to an alternaitve genome\n";	
		print $debug_print "-> Keep read\n";
	}
	
}


sub sum_matches {
	my $cigar = shift;
	my @matched;
	push @matched, $cigar =~ /.*?(\d+)M/g;

	my $match_score;
	$match_score += $_ foreach @matched;
	
	return ($match_score, @matched);
}

sub usage {
	say "******** cigar_compare.pl ********";
    say "Usage: $0 [options]";
	say "  --bam = specify bam file processed with 'tag_reads.py'";
	say "  --softclip = filter reads where both ends are softclipped > 20 and have a match < 30";
	say "  --test  = perform run on test dataset 'files/test.bam'";
	say "  --debug = run in debug mode";
	say "  --help = print this help message and exit\n";
	say "Nick Riddiford 2017";
}