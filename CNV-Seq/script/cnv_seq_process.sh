#!/bin/bash

Rscript='/data/kdi_prod/project_result/948/01.00/Analysis/Analysis/CNV-Seq/script/cnv-seq.R'


for file in $@;
do
  stem=$(basename "${file}" )
  id=$(echo $stem | cut -d'.' -f 1)
  
  Rscript $Rscript $file $id --save
done

