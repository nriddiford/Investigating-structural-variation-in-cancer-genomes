#!/usr/bin/perl
use strict;
use warnings;

use Data::Dumper;


my %seqs;
my $head;

my %lengths;

my $longest_contig = 0;
while(<>){
  chomp;

  if (/^>/){
    $head = $_;
  }

  my $seq = $_ if /^[ATGC]/;


  if ($head and $seq){
    $seqs{$head} .= $seq ;

    if ( (not $longest_contig) or ($longest_contig < length $seqs{$head}) ){
      $longest_contig = length $seqs{$head};
    }

  }

}

print $longest_contig;


__DATA__
>Contig1
AAGTGGTAGTATCAAATGTAATGTACTCTTTCATCAATAGATTTAGGGCCGATAGAATGG
TAAAAGTTGTAGTTTCAATGGAATAATGTGCAATATACGGCACTTGACCAT
>Contig2
GCAATATACGGCACTTGACCATTGACATAGCTATCAAGAATTGAATTCACATTCAATTTG
AGGCATGTCTCGAGAAATTGATCGGCCAAAACTTGAATTGCTGGCGTTGTGTTGGCTGCA
CATTTGATACCGCCCTCGGGATTGATAATG
>Contig3
TTTGATACCGCCCTCGGGATTGATAATGTACAATTCACCGCGACATGCCCAAATGTACCA
ATCGAACGGTGAGACGACAATTAAACGATTGCATACCGTACTTGTGGATTTGTCGTAGAT
TTTGAAAATTGTAGGTGTAAAAGCTGATTGTTGGGC
