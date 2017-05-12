#!/bin/bash

cn2bed=/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/CNV-Seq/script/CN2bed.py
bed_sort=/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/CNV-Seq/script/bed_sort.py
cnv2gff=/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/CNV-Seq/script/cnv2gff.pl

for file in $@; do
  stem=$(basename "${file}" )
  id=$(echo $stem | cut -d'_' -f 1)

  echo "Converting CNVs to bedpe for $file"

  python $cn2bed -c $file -o ${id}.bedpe
  python $bed_sort -b ${id}.bedpe > ${id}_cnvs.bedpe
  
  rm ${id}.bedpe

  echo "Creating gff track"
  perl $cnv2gff $file

done
