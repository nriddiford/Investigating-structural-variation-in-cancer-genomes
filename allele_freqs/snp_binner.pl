#!/usr/bin/perl
use strict;
use warnings;

use feature qw/ say /;

use Data::Printer;

use File::Basename;

use FindBin;
use FindBin '$Script';

die usage() unless @ARGV;

my $file = $ARGV[0];

open my $in, '<', $file or die $!;
 
my (@file_parts) = split(/\./, basename($file));

my ($id) = $file_parts[0];

open my $out, '>', $id . '.allele_freqs.txt' or die $!;

print $out join("\t", "chrom", "position", "freq", "type", "variants_included") . "\n";

my %karyotype = qw / 2L 23513712 2R 25286936 3L 28110227 3R 32079331 4 1348131 X 23542271 Y 3667352 /;

my %snvs;
while(<$in>){
	chomp;
	next if /chrom/;
	my ($chrom, $pos, $n_freq, $t_freq) = (split)[0,1,6,10];
	($n_freq) =~ s/\%//;
	($t_freq) =~ s/\%//;
	next unless $karyotype{$chrom};
	
	$snvs{$chrom}{$pos} = [$n_freq, $t_freq];
	
}	

my $window_size = 100000;

my %freq_by_window;

for my $chrom (sort keys %karyotype){
	my ($n_freq, $t_freq);
	
	my $window_count = 0;
	
	my ($w_start, $w_stop) = (0,0);
	
	for (my $i = 0; $i <= $karyotype{$chrom}; $i += $window_size) {

		my $w_start = $i;
		my $w_stop = $w_start + ($window_size - 1);
				
		my ($cum_n_freq, $cum_t_freq, $row_count) = (0,0,0);
		
		for my $pos (keys %{$snvs{$chrom}} ){
			($n_freq, $t_freq) = @{$snvs{$chrom}{$pos}};

			if ( $pos >= $w_start and $pos < $w_stop){
				$row_count++;
				$cum_n_freq += $n_freq;
				$cum_t_freq += $t_freq;
			}

		}

		if ( $row_count >= 20 ){
			my $location = ( $w_start + ($window_size / 2 ));
			
			# print join ("\t", $chrom, $w_start . "-" . $w_stop, $n_freq, $t_freq, "row_count:", $row_count) . "\n";

			my ($av_n_freq) = ($cum_n_freq + 0.01)/$row_count;
			my ($av_t_freq) = ($cum_t_freq + 0.01)/$row_count;

			print $out join ( "\t", $chrom, $location, $av_n_freq, "normal_freq", $row_count ) . "\n";
	 		print $out join ( "\t", $chrom, $location, $av_t_freq, "tum_freq", $row_count ) . "\n";
		}
		
	}
	
}


sub usage {
	say "********** $Script ***********";
	say "perl $Script <allele.freqs.txt>";
}