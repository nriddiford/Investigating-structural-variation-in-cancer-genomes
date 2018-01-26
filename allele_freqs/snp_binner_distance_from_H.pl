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

open my $out, '>', $id . '.dist_from_h.txt' or die $!;

print $out join("\t", "chrom", "position", "window", "dist_ratio", "count") . "\n";

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
				
		my ($cum_n_dist, $cum_t_dist, $row_count) = (0,0,0);
		
		for my $pos (keys %{$snvs{$chrom}} ){
			($n_freq, $t_freq) = @{$snvs{$chrom}{$pos}};
			my ($normal_distance_from_h, $tumour_distance_from_h) = (0,0);
			
			if ( $pos >= $w_start and $pos < $w_stop){
				$row_count++;
				
				$normal_distance_from_h = abs(50 - $n_freq);
				$tumour_distance_from_h = abs(50 - $t_freq);
			
				$cum_n_dist += $normal_distance_from_h;
				$cum_t_dist += $tumour_distance_from_h;
				
				# print join ("\t", $chrom, $w_start . "-" . $w_stop, $n_freq, $t_freq, "n_dist: ", $normal_distance_from_h, "t_dist: ", $tumour_distance_from_h,  "row_count:", $row_count) . "\n";
							
			}

		}
		
		if ( $row_count >= 20 ){
			my $location = ( $w_start + ( $window_size / 2 ));
			
			my $av_n_dist = $cum_n_dist / $row_count;
 			my $av_t_dist = $cum_t_dist / $row_count;
			
			my $dist_ratio = $av_t_dist/$av_n_dist;
			
			print $out join ("\t", $chrom, $location, $w_start . "-" . $w_stop, $dist_ratio, $row_count) . "\n";
			
			# print $out join ("\t", $chrom, $location, $w_start . "-" . $w_stop, $av_n_dist, "normal", $row_count) . "\n";
			# print $out join ("\t", $chrom, $location, $w_start . "-" . $w_stop, $av_t_dist, "tumour", $row_count) . "\n";

		}
		
	}
	
}


sub usage {
	say "********** $Script ***********";
	say "perl $Script <allele.freqs.txt>";
}