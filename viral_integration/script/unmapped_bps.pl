#!/usr/bin/perl

use strict;
use warnings;

use feature qw/ say /;

if ( $#ARGV != 1 ){
  die "Usage clusters.bed unmapped_reads.bam";
}

my $file = $ARGV[0];

open my $in, '<', $file or die "$! $file";

my $cluster=0;

my $bam = $ARGV[1];

while(<$in>){
  chomp;

  $cluster++;

  my ($chrom, $cluster_start, $cluster_end) = split;

  my $start = ($cluster_start - 100);
  my $end = ($cluster_end + 100);

  say "Extracting unmapped reads in cluster $cluster from $bam";
  say "  Cluster $cluster region: $chrom:$start-$end";
  
  system("samtools view -f 4 -F 264 $bam $chrom:$start-$end -o cluster\_$cluster\_$chrom\_$start\_$end.bam");

  say "Getting fasta for cluster $cluster";

  system("python script/bam2fa.py cluster\_$cluster\_$chrom\_$start\_$end.bam");
  system("rm cluster\_$cluster\_$chrom\_$start\_$end.bam");

  say "Done for cluster $cluster";

}

__DATA__
2L      1244608 1244823
2L      1301135 1301236
2L      1303777 1303878
2L      1655480 1655581
2L      1789278 1789379
2L      1964091 1964192
2L      2729641 2729740
2L      5071650 5071751
2L      5870803 5870967
2L      6203119 6203163
2L      6403041 6403138
2L      6801377 6801478
2L      6801512 6801613
2L      7397954 7398110
2L      7714093 7714194
