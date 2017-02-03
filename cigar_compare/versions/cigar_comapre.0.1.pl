#!/usr/bin/perl

use strict;
use warnings;

use feature qw/ say /;

use Getopt::Long qw/ GetOptions /;

my $debug;
my $help;
my $test;
my $bam_file = 'files/test.bam';

GetOptions( 'bam=s'				=>			\$bam_file,
			'debug'				=>			\$debug,
			'test'				=>			\$test, 
			'help'				=>			\$help
	  	  )	or die usage();
		  
if ($help)  { exit usage() } 

my $in;

if ($test) {
	open $in , '<', $bam_file or die $!;
}
else {
	open $in, "samtools view -h $bam_file |" or die $!;
}

while(<$in>){
	chomp;
	
	if (/^@/){
		say unless $debug;
		next;
	}

	my ($read) = (split)[0];
	my ($droso_cigar) = (split)[5];
	next unless $droso_cigar =~ /M/;
	my ($droso_matched) = $droso_cigar =~ /.*?(\d+)M/;

	if (/CIGAR:/){
	
		my $alt_cigar_group;
		my ($alt_cigar, $alt_matched);
		my @split = split;
		
		my $alt_mapping_flag = 0;
		for (@split){
			if (/^[CE]D/){
				$alt_cigar_group = $_;
				($alt_cigar, $alt_matched) = $alt_cigar_group =~ /CIGAR:(.*?(\d+)M)/;
				$alt_mapping_flag = 1;
			}
		}
		
		if ($alt_mapping_flag){
			
			if ($debug){
				say "_______________________________";
				say $read;
	 			print "Droso:    $droso_cigar    $droso_matched matched\n";
	 			print "Alt:    $alt_cigar    $alt_matched matched\n";
				($droso_matched > $alt_matched) ? say "-> Keep: $droso_matched > $alt_matched" : say "-> Throw: $droso_matched < $alt_matched";
			}
			
			else {
				($droso_matched > $alt_matched) ? say : next;
			}
			
		}
	
	 }
	 else { say unless $debug }
}


sub usage {
	say "******** cigar_compare.pl ********";
    say "Usage: $0 [options]";
	say "  --bam = specify bam file processed with 'tag_reads.py'";
	say "  --debug = run in debug mode";
	say "  --test  = run in test mode";
	say "  --help = print this help message and exit\n";
	say "Nick Riddiford 2017";
}