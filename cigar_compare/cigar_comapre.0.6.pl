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

my ($debug_print, $sc_filter_print, $not_filtered_bed);
if ($debug){
	open $debug_print, '>', $id . ".debug.log" or die $!; 
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
	
	# For bed output
	my ($chrom, $start, $stop, $seq) = (split)[2,3,7,9];
	my $read_length = length($seq);
	my $bed_stop = $start + $read_length;
	
	
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
		my $non_filter_alt_cigar_group;
		my ($alt_cigar, $alt_matched);

		my $alt_mapping_flag = 0;
		my $non_filter_set = 0;
				
		for (split){
			if (/^[CEG]D/){
				($alt_cigar_group) = $_;
				$alt_mapping_flag = 1;
			}
			elsif (/^[A]D/){
				$non_filter_alt_cigar_group = $_;
				$non_filter_set = 1;
			}
		}
		
		my ($non_filt_alt_cigar, $non_filt_alt_matched);
		
		
		if ($non_filter_set){

			($non_filt_alt_cigar, $non_filt_alt_matched) = $non_filter_alt_cigar_group =~ /CIGAR:(.*?(\d+)M.*?),/g;
		
			my @non_filt_matched;
			push @non_filt_matched, $non_filt_alt_cigar =~ /.*?(\d+)M/g;
			say unless $debug;
			print $debug_print "Alt non-filter:    $non_filt_alt_cigar " . join(', ', @non_filt_matched) . " matched\n" if $debug;
			print $debug_print "-> Keep: $droso_match_score [Alternative alignment in non-filter set]\n" if $debug;
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
		elsif ($droso_match_score < $alt_match_score){
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
	
	elsif ($soft_clip and $droso_cigar =~ /\d+S\d+M\d+S/){	
		my ($sc_1, $m, $sc_2) = $droso_cigar =~ /(\d+)S(\d+)M(\d+)S/;
			
			if ($sc_1 > 15 and $m < 30 and $sc_2 > 15){
				print $sc_filter_print "$chrom\t$start\t$bed_stop\n" if $debug;
				print $debug_print "-> Throw: $droso_cigar [exceeded double-ended soft-clipping threshold]\n" if $debug;
				next;
			}
			
	}
	
	else {
		print $debug_print "-> Keep: $droso_match_score [no alt mapping]\n" if $debug;
		say unless $debug;
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
	say "  --softclip = filter reads where both ends are softclipped > 20 and have a match < 30";
	say "  --test  = perform run on test dataset 'files/test.bam'";
	say "  --debug = run in debug mode";
	say "  --help = print this help message and exit\n";
	say "Nick Riddiford 2017";
}