#!/usr/bin/perl

use strict;
use warnings;

use feature qw/ say /;

use Getopt::Long qw/ GetOptions /;

my $debug;
my $help;
my $test;
my $bam_file = $ARGV[0];

GetOptions( 'bam=s'				=>			\$bam_file,
			'debug'				=>			\$debug,
			'test'				=>			\$test, 
			'help'				=>			\$help
	  	  )	or die usage();
		  
if ($help)  { exit usage() } 

open my $in, '<', $bam_file or die $!;

while(<$in>){
	chomp;
	if (/^@/){
		say unless $debug;
		next;
	}
	
	my $read = (split)[0];
	my ($droso_cigar) = (split)[5];
	my ($alt_cigar_group) = (split)[11];
	my ($droso_matched) = $droso_cigar =~ /.*?(\d+)M/;	
	my ($alt_cigar, $alt_matched) = $alt_cigar_group =~ /CIGAR:(.*?(\d+)M)/;
	
	# my ($alt_mate_cigar_group) = (split)[12];
	# my ($alt_mate_cigar, $alt_mate_matched) = $alt_mate_cigar_group =~ /CIGAR:(.*?(\d+)M)/;
	
	if ($debug){
		say "_______________________________";
		say $read;
 		print "Droso:    $droso_cigar    $droso_matched matched\n";
 		print "Alt:    $alt_cigar    $alt_matched matched\n";
		($droso_matched > $alt_matched) ? say "-> Keep: $droso_matched > $alt_matched" : say "-> Throw: $droso_matched < $alt_matched";
	}
	
	unless ($debug){
		($droso_matched > $alt_matched) ? say : next;
	}
	
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