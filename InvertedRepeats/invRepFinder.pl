#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use File::Basename;
use Data::Dumper;
use feature qw/ say /;

my $genome_file = '/Users/Nick_curie/Documents/Curie/Data/Genomes/dmel_6.12.fa';# work
my $genome_file = '/data/kdi_prod/project_result/948/01.00/Analysis/Genomes/Dmel_6/dmel_6.12.fa'; # cluster
#my $genome_file = '/Users/Nick/Documents/Curie/Data/Genomes/Dmel_v6.12/Dmel_6.12.fasta'; # home

my %chroms = qw / 2L 23513712 2R 25286936 3L 28110227 3R 32079331 4 1348131 X 23542271 Y 3667352 /;

my $genome_ref = get_genome($genome_file);
my %genome = %{$genome_ref};

sub get_genome {
  my $genome_file = shift;
  my $seqio = Bio::SeqIO->new('-file' => "$genome_file", '-format' => 'Fasta');

  my $genome_length;

  say "Reading in genome: $genome_file";

  while(my $seq = $seqio->next_seq) {
    my $nucs = $seq->seq;
    my $chr = $seq->id;
    next unless $chroms{$chr};
    $genome{$chr} = $nucs;
  }
  return(\%genome);
}


#print substr($dna, 106, 10) . "\n";


open my $sir_out, '>', 'sir_positions.txt';
open my $cru_out, '>', 'cruciform_positions.txt';
open my $lir_out, '>', 'lit_positions.txt';
open my $bed_out, '>', 'repeats.bed';

# From Woczic et al 2012
# Two tracts ≥10 nt with same sequence on the complementary strands, separated by up to half tract length but no more than 20 nt

my $min = 10;
my $max = 100;
my $min_spacer = 0;
my $max_spacer = 2000;

for my $chromosome (sort keys %genome ){
  say "Searching for inverted repeats on chromosome: $chromosome";
  say "Repeats between $min and $max bp and spacer between $min_spacer and $max_spacer";

  my $dna_length= length($genome{$chromosome});
  my ($length, $pos, $region, $search_window, $inverted_region);

  for ($pos = 0 ; $pos <= $dna_length-$max-$max-$max_spacer ; $pos++) {
    for ($length = $min ; $length <= $max ; $length++) {
      $region = substr($genome{$chromosome}, $pos, $length);
      $search_window = substr($genome{$chromosome}, $pos-1, 3000);
      $inverted_region = reverse $region;
      $inverted_region =~ tr/[ACGT]/[TGCA]/;

        if ($search_window =~ /^[ACGT](($region)([ACGT]{$min_spacer,$max_spacer})($inverted_region))/) {

          # Capture groups:
          # $1 = entire sequence
          # $2 = region seq
          # $3 = spacer seq
          # $4 = inv_region seq

          my $igv_start = $pos + 1;
          my $igv_end = $igv_start + length($1) - 1;
          my $zero_based_start = $pos;
          my $end = $zero_based_start + length($1);

          my $size = length($2);
          my $e2_start = length($1) - $size;


          if ( (length($3) < length($region)/2) and (length($3) < 20) ){
            print "Cruciform repeat starts at $chromosome:$igv_start => $2***$3***$4\n$1\n";
            print "$chromosome:$igv_start-$igv_end\n\n";
            print $cru_out "$chromosome\t$zero_based_start\n";
            print $cru_out "$chromosome\t$end\n";
            # Red track
            print $bed_out join(" ", $chromosome, $zero_based_start, $end, "Cruciform", '10', '+', $zero_based_start, $end, '255,0,0', '2', "$size,$size", "0,$e2_start") . "\n";
          }
          elsif (length($1) > 500 and length($region) >= 30){
            print "Long inverted repeat starts at $chromosome:$igv_start => $2***$3***$4\n$1\n";
            print "$chromosome:$igv_start-$igv_end\n\n";
            print $lir_out "$chromosome\t$zero_based_start\n";
            print $lir_out "$chromosome\t$end\n";
            # Blue track
            print $bed_out join(" ", $chromosome, $zero_based_start, $end, "LIR", '10', '+', $zero_based_start, $end, '0,0,255', '2', "$size,$size", "0,$e2_start") . "\n";

          }

          else{
          print "Short inverted repeat starts at $chromosome:$igv_start => $2***$3***$4\n$1\n";
          print "$chromosome:$igv_start-$igv_end\n\n";
          print $sir_out "$chromosome\t$zero_based_start\n";
          print $sir_out "$chromosome\t$end\n";
          print $bed_out join(" ", $chromosome, $zero_based_start, $end, "SIR", '10', '+', $zero_based_start, $end, '0', '2', "$size,$size", "0,$e2_start") . "\n";
          # Black track
          say "Bed file coordinates:";
          print join(" ", $chromosome, $zero_based_start, $end, "SIR", '10', '+', $zero_based_start, $end, '1', '2', "$size,$size", "0,$e2_start") . "\n";

        }

      }
    }
  }
}
