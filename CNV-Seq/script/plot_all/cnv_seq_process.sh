#!/bin/bash

Rscript='/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/CNV-Seq/script/cnv-seq.R'
#chrom='2L'

for chrom in 2L 2R 3L 3R 4 X Y; do
  for file in $@; do
    
    stem=$(basename "${file}" )
    id=$(echo $stem | cut -d'.' -f 1)
    
    echo "Plotting $chrom for $id"

    Rscript $Rscript $file $id $chrom --save
  done

  rm Rplots.pdf

  if [ ! -d plots ]; then
    mkdir -p plots;
  fi

  mv *pdf plots/

done
