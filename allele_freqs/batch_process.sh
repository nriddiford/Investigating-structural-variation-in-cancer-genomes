#!/bin/bash

Rscript='line_plot.r'

if [ ! -d plots ]; then
  mkdir -p plots;
fi

for file in "$@"; do

  stem=$(basename "${file}" )
  id=$(echo $stem | cut -d'.' -f 1)

  echo "Plotting all chroms for $id"
  
  ./$Rscript $file
  
done
 
  # mv *pdf plots/
  #
  # rm Rplots.pdf