#!/bin/bash

Rscript='cnv-seq.R'
chrom='X'

for file in $@; do
        stem=$(basename "${file}" )
        id=$(echo $stem | cut -d'.' -f 1)
        Rscript $Rscript $file $id $chrom --save
done

rm Rplots.pdf

if [ ! -d plots ]; then
    mkdir -p plots;
fi

mv *pdf plots/