# Pipeline

# Table of Contents
* [SVDetect](#svdetect)
  * [About the tool](#about-the-tool)
  * [Protocol](#protocol)
* [Trimming and mapping](#trimming-and-mapping)
  * [FastQC](#run-fastqc-to-visually-inspect-all-sequencing-results)
  * [Trimmomatic](#run-trimmomatic-to-quality-filter-reads)
  * [Bowtie index](#build-bowtie-index)
  * [Tophat](#run-tophat2-on-each-set-of-paired-reads)
* [Transcript assembly and differential expression analysis](#transcript-assembly-and-differential-expression-analysis)
  * [Cufflinks](#assemble-transcripts-and-calcuulate-abundance-estimation-with-cufflinks)
  * [Cuffmerge](#merge-assemblies-using-cuffmerge)
  * [Cuffdiff](#differential-expression-using-cuffdiff)
* [Estimating variance between biological replicates](#estimating-variance-between-biological-replicates)
* [Annotating transcript models](#annotating-transcript-models)
  * [Build transcript models](#build-transcript-models)
  * [Add differential expression information](#add-differential-expression-information)
  * [Annotate](#annotate)
* [Gene Enrichment Analysis](#gene-enrichment-analysis)
* [Find co-differentially expressed genes](#find-co-differentially-expressed-genes)
  * [For single conditions](#for-single-conditions)
  * [For merged conditions](#for-merged-conditions)
* [Gene Ontology](#gene-ontology-analysis-on-discrete-gene-sets)
  * [DAVID](#david)



# SVDetect

### About the tool
[SVDetect](http://bioinformatics.oxfordjournals.org/content/26/15/1895.full.pdf) is a tool to identify genomic structural variations from paired-end and mate-pair sequencing data
Applying both sliding-window and clustering strategies, it uses anomalously mapped read pairs to localise genomic rearrangements and classify them according to their type e.g.:
* Large insertions and deletions
* Inversions
* Duplications
* Balanced and unbalanced inter-chromosomal translocations


### Protocol

