#!/usr/bin/perl

use strict;
use warnings;

use feature qw/ say /;

use Data::Dumper;

use Getopt::Long qw/ GetOptions /;


my $individual_file;

GetOptions( 'file=s'					=>			\$individual_file
		   ) or die usage();
		 
	
my @files;
	 
if ($individual_file){
	push @files, $individual_file;
}

else {
	@files = qw/ lumpy_svs.txt lumpy_bnds.txt cnv-seq_calls.txt delly_bnds.txt delly_svs.txt freec_calls.txt /;
}

my %karyotype = qw / 2L 23513712 2R 25286936 3L 28110227 3R 32079331 4 1348131 X 23542271 Y 3667352 /;

my $window_size = 1000000;
$window_size = ($window_size - 1);
my (@chroms, @bp1, @bp2s);
my %SVs;
my $file_count = 0;

my $sv_id = 0;
for (@files){
	open my $file, '<', $_ or die $!;
	say "Reading '$_'";
	process($file, $sv_id);

}

my %SVs_per_window;
	
for my $chrom (keys %karyotype){
	my ($bp1, $bp2);
	
	my $j = 0;
	my $window_count = 0;
	
	for (my $i = 0; $i <= $karyotype{$chrom}; $i += $window_size) {
		$window_count++;
		$j = $i + $window_size;
		my $sv_count = 0;
		for my $var (keys %{$SVs{$chrom}} ){
			($bp1, $bp2) = @{$SVs{$chrom}{$var}};
			
			if ($bp1 <= $j and $bp1 >= $i){
				$sv_count++;
				$SVs_per_window{$chrom}{$window_count} = [$sv_count, $i, $j] ;
			}
			
			next if $bp1 == $bp2;
			
			if ($bp2 <= $j and $bp2 >= $i){
				$sv_count++;
				$SVs_per_window{$chrom}{$window_count} = [$sv_count, $i, $j] ;
			}
			
		}
		
		if ($sv_count == 0){
			$SVs_per_window{$chrom}{$window_count} = [0, $i, $j];
		}
		
	$i++;
	}
}

open my $out, '>', 'all_structural_variants.txt' or die $! if not $individual_file;
open my $g4s, '>', 'G4s.txt' or die $! if $individual_file;


my $total_count = 0;
for my $chr (sort keys %SVs_per_window){
	for my $w (sort { $a <=> $b } keys $SVs_per_window{$chr}){
		my ($count, $start, $stop) = @{$SVs_per_window{$chr}{$w}};
		print $out join("\t", $chr, $start, $stop, $count) . "\n" if not $individual_file;
		print $g4s join("\t", $chr, $start, $stop, $count) . "\n" if $individual_file;
		$total_count += $count;
	}
}

if ($individual_file){
	say "Writing all G4 calls to 'G4s.txt' [G4 quad heatmap file]";
	
	say "-------------------------------------------";
	say " $total_count G4 sequences detected ";
	say "-------------------------------------------";
}

else{
	say "Writing all SV calls to 'all_structural_variants.txt' [SV heatmap file]";

	say "-------------------------------------------";
	say " $total_count motifs detected ";
	say "-------------------------------------------";
}


sub process { 
	my ($file, $count) = @_;
	
	while(<$file>){
		chomp;
		my ($chrom, $bp1, $bp2) = (split)[0..2];
		$sv_id++;
		$SVs{$chrom}{$sv_id} = [$bp1, $bp2];
		
	}
	return %SVs;
}

