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

print $out join("\t", "chrom", "window", "position", "freq", "type", "change_dir", "variants_included") . "\n";

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
				
		my ($cum_n_freq, $cum_t_freq, $tum_inc_count, $tum_dec_count, $row_count) = (0,0,0,0,0);
		
		my ($cum_n_high, $cum_t_high, $cum_n_low, $cum_t_low) = (0,0,0,0);
		
		for my $pos (keys %{$snvs{$chrom}} ){
			($n_freq, $t_freq) = @{$snvs{$chrom}{$pos}};

			if ( $pos >= $w_start and $pos <= $w_stop){
				$row_count++;
				
				if ($n_freq < $t_freq){
					$cum_n_low += $n_freq;
					$cum_t_high += $t_freq;
					$tum_inc_count++;
				}
				
				elsif ($n_freq > $t_freq) {
					$cum_n_high += $n_freq;
					$cum_t_low += $t_freq;
					$tum_dec_count++;
				}
				
			}

		}

		if ( $row_count >= 20 ){
			my $position = ( $w_start + ($window_size / 2 ));
			
			if ($tum_dec_count){
				
				# print join ( "\t", $chrom, $w_start . "-" . $w_stop, $cum_n_high, "normal", "n>t", $tum_dec_count ) . "\n";
				
				my ($av_n_high_freq) = ($cum_n_high + 0.01) / $tum_dec_count;
				
				print $out join ( "\t", $chrom, $w_start . "-" . $w_stop, $position, $av_n_high_freq, "normal", "n>t", $tum_dec_count ) . "\n";
				
				# print join ( "\t", $chrom, $w_start . "-" . $w_stop, $cum_t_low,  "tumour", "n>t", $tum_dec_count ) . "\n";
				
				my ($av_t_low_freq) =  ($cum_t_low  + 0.01) / $tum_dec_count;
				
				print $out join ( "\t", $chrom, $w_start . "-" . $w_stop, $position, $av_t_low_freq,  "tumour", "n>t", $tum_dec_count ) . "\n";
					
			}
			
			if ($tum_inc_count){
				
				# print join ( "\t", $chrom, $w_start . "-" . $w_stop, $cum_n_low,  "normal", "t>n", $tum_inc_count ) . "\n";
				
				my ($av_n_low_freq) = ($cum_n_low + 0.01)/$tum_inc_count;
				
				print $out join ( "\t", $chrom, $w_start . "-" . $w_stop, $position, $av_n_low_freq,  "normal", "t>n", $tum_inc_count ) . "\n";
				
				my ($av_t_high_freq) = ($cum_t_high + 0.01)/$tum_inc_count;
								
				# print join ( "\t", $chrom, $w_start . "-" . $w_stop, $cum_t_high, "tumour", "t>n", $tum_inc_count ) . "\n";
				
				print $out join ( "\t", $chrom, $w_start . "-" . $w_stop, $position, $av_t_high_freq, "tumour", "t>n", $tum_inc_count ) . "\n";
				
			}

		
		}
		
	}
	
}


sub usage {
	say "********** $Script ***********";
	say "perl $Script <allele.freqs.txt>";
}