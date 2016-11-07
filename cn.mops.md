# cn.mops protocol

# Table of Contents
* [About the tool](#about-the-tool)
* [Protocol](#protocol)


## About the tool

## Protocol

```{R}
library(cn.mops)
setwd("/data/kdi_prod/project_result/948/01.00/Analysis/Trimmo_out/Bwa_test/Bwa_full")

BAMFiles <- list.files(system.file("extdata", package="cn.mops"),pattern="sorted.bam$",
+                 full.names=TRUE)
bamDataRanges <- getReadCountsFromBAM(BAMFiles,
+                 sampleNames=paste("Sample",1:3),mode="paired")



BAMFiles <- list.files(pattern="sorted.bam$")
bamDataRanges <- getReadCountsFromBAM(BAMFiles,mode="paired")
res <- cn.mops(bamDataRanges)
```
