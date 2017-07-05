#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l walltime=02:00:00
#PBS -l mem=1GB
#PBS -N clusters
#PBS -o /data/kdi_prod/project_result/948/01.00/Analysis/Analysis/Characterised_SVs/viral_int
#PBS -e /data/kdi_prod/project_result/948/01.00/Analysis/Analysis/Characterised_SVs/viral_int

# HUM
# A512
# A558
# A572
# A370

group=HUM
bam_files=/data/kdi_prod/project_result/948/01.00/Analysis/Bwa/${group}
out_dir=/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/clusters
script_bin=/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/clusters/script

if [ ! -d $out_dir/${group} ]
then
    mkdir -p $out_dir/${group}
fi

if [ ! -d $out_dir/${group}/bam_out ]
then
    mkdir -p $out_dir/${group}/bam_out
fi

if [ ! -d $out_dir/${group}/fasta ]
then
    mkdir -p $out_dir/${group}/fasta
fi

if [ ! -d $out_dir/${group}/bed ]
then
    mkdir -p $out_dir/${group}/bed
fi

for file in $(ls -1 "$bam_files"/*.tagged.SC.RG.bam | sort -V )
do # sort -V doesn't work for local...

  stem=$(basename "${file}" )
  id=$(echo $stem | cut -d '.' -f 1)

  if [ ! -f $out_dir/${group}/bam_out/${id}.mate_not_mapped.bam ]
  then
    echo "Extracting reads with unmapped mates from $stem"
    echo "samtools view -q 30 -f 8 -F 1284 $file -o ${out_dir}/${group}/bam_out/${id}.mate_not_mapped.bam"
    samtools view -q 30 -f 8 -F 1284 $file -o ${out_dir}/${group}/bam_out/${id}.mate_not_mapped.bam
  fi

  if [ ! -f $out_dir/${group}/bam_out/${id}.read_unmapped_mate_mapped.bam ]
  then
    echo "Extracting unmapped reads whose mates are mapped and primary from $stem"
    echo " samtools view -f 4 -F 264 $file -o ${out_dir}/${group}/bam_out/${id}.read_unmapped_mate_mapped.bam"
    samtools view -f 4 -F 264 $file -o ${out_dir}/${group}/bam_out/${id}.read_unmapped_mate_mapped.bam
  fi

done

set -- $(ls -1 $out_dir/${group}/bam_out/*.mate_not_mapped.bam | sort -V ) # sort -V doesn't work for local...

[[ -e $1 || -L $1 ]] || { echo "No mate_not_mapped.bam files found in current dir" >&2; exit 1; }

export PATH=$HOME/miniconda2/bin:$PATH
source activate unmapped

while (( $# > 1 )); do
  tum_stem=$(basename "${1}" )
  tum_id=$(echo $tum_stem | cut -d '.' -f 1)

  con_stem=$(basename "${2}" )
  con_id=$(echo $con_stem | cut -d '.' -f 1)

  echo "Writing bed file for merged clusters of 5 or more reads for TUM: $tum_id"
  bedtools bamtobed -i ${out_dir}/${group}/bam_out/${tum_id}.mate_not_mapped.bam | bedtools cluster -s -d 100 > ${out_dir}/${group}/bed/${tum_id}.mate_not_paired_clusters.bed

  awk 'NR==FNR{a[$7]++;next}a[$7]>4' ${out_dir}/${group}/bed/${tum_id}.mate_not_paired_clusters.bed ${out_dir}/${group}/bed/${tum_id}.mate_not_paired_clusters.bed | bedtools sort | bedtools merge -d 100 > ${out_dir}/${group}/bed/${tum_id}.clusters.bed

  echo "Writing bed file for merged reads for NORM: $con_id"
  bedtools bamtobed -i ${out_dir}/${group}/bam_out/${con_id}.mate_not_mapped.bam | bedtools merge -d 100 > ${out_dir}/${group}/bed/${con_id}.mate_not_paired_merged.bed

  echo "Subtracting clusters from ${out_dir}/${group}/bed/${tum_id}.clusters.bed containing reads from NORM: $con_id"
  bedtools subtract -A -a ${out_dir}/${group}/bed/${tum_id}.clusters.bed -b ${out_dir}/${group}/bed/${con_id}.mate_not_paired_merged.bed > $out_dir/${group}/${tum_id}.clusters_sub.bed

  echo "Indexing ${out_dir}/${group}/bam_out/${tum_id}.read_unmapped_mate_mapped.bam"

  samtools index ${out_dir}/${group}/bam_out/${tum_id}.read_unmapped_mate_mapped.bam

  echo "Getting fasta for unmapped reads surrounding clusters in ${tum_id}.clusters_sub.bed"

  echo "perl $script_bin/unmapped_bps.pl $out_dir/${group}/${tum_id}.clusters_sub.bed $out_dir/${group}/bam_out/${tum_id}.read_unmapped_mate_mapped.bam"

  perl $script_bin/unmapped_bps.pl $out_dir/${group}/${tum_id}.clusters_sub.bed $out_dir/${group}/bam_out/${tum_id}.read_unmapped_mate_mapped.bam

  rename "s/cluster_/${tum_id}.cluster_/g" $out_dir/*.fa


  if [ ! -d $out_dir/${group}/fasta/${tum_id} ]
  then
    mkdir -p $out_dir/${group}/fasta/${tum_id}
  fi

  mv $out_dir/*.fa $out_dir/${group}/fasta/${tum_id}

  for fasta_file in $out_dir/${group}/fasta/${tum_id}/${tum_id}*.fa
  do
    cap3 $fasta_file -o 16 -p 75 -z 2 -s 500

    file_name=$( basename "${fasta_file}" )
    cluster_id=$(echo $file_name | cut -d '.' -f 1,2)

      if [ $(perl $script_bin/get_contigs.pl $out_dir/${group}/fasta/${tum_id}/${cluster_id}.fa.cap.contigs ) -lt 130 ]
      then
         rm $out_dir/${group}/fasta/${tum_id}/${cluster_id}.*

      elif [ ! -s "$out_dir/${group}/fasta/${tum_id}/${cluster_id}.fa.cap.contigs" ]
      then
        rm $out_dir/${group}/fasta/${tum_id}/${cluster_id}.*
      fi

  done

  shift 2 || break

done
