#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use File::Basename;

use Data::Dumper;

use Getopt::Long qw/ GetOptions /;

my $vcf_file;
my $pon_file;
my $help;

GetOptions( 'vcf=s'                  =>    \$vcf_file,
            'pannel-of-normals=s'    =>    \$pon_file,
            'help'                   =>    \$help
          ) or die usage();

if ($help) { exit usage() }

# my $keep_file = shift;
# my $filt_file = shift;


open my $somatic_calls, '<', $vcf_file;
open my $germ_calls, "gunzip -c $pon_file |";

my %germlineSNVs;
while(<$germ_calls>){
  chomp;
  next if /^#/;
  my ($chrom, $pos, $alt) = (split)[0,1,4];
  $germlineSNVs{$chrom}{$pos}{$alt} = 1;
}

my ($name, $extention) = split(/\.([^.]+)$/, basename($vcf_file), 2);

my ($id) = (split(/_/, $name))[0];


open my $pon_filtered, '>', $id . '_filt_somatic.vcf';

my $germline_filtered = 0;
my $remaining_calls = 0;
while(<$somatic_calls>){
  chomp;
  if (/^#/){
    print $pon_filtered "$_\n";
    next;
  }
  my ($chrom, $pos, $alt) = (split)[0,1,4];

  if ($germlineSNVs{$chrom}{$pos}{$alt}){
    $germline_filtered++;
    next;
  }
  $remaining_calls++;
  print $pon_filtered "$_\n";
}

print "$germline_filtered germline calls filtered\n$remaining_calls unfiltered calls\n";



sub usage {
  print
"
usage: $0 [-h] [-v VCF] [-p PON]

varscanfilt
author: Nick Riddiford (nick.riddiford\@curie.fr)
version: v1.0
description: Filter Varscan VCF file against a pannel of normals created by Mutect2 (or similar)

arguments:
  -h, --help            show this help message and exit
  -v VCF, --vcf
                        VCF file from Varscan2[required]
  -p PON --pannel-of-normals
                        PON file (e.g. created from Mutect2)
"
}
