# SVDetect protocol

# Table of Contents
* [About the tool](#about-the-tool)
* [Protocol](#protocol)
* [Output](#output)

## About the tool
[SVDetect](http://bioinformatics.oxfordjournals.org/content/26/15/1895.full.pdf) is a tool to identify genomic structural variations from paired-end and mate-pair sequencing data.
Applying both sliding-window and clustering strategies, it uses anomalously mapped read pairs to localise genomic rearrangements and classify them according to their type e.g.:
* Large insertions and deletions
* Inversions
* Duplications
* Balanced and unbalanced inter-chromosomal translocations

Starting from a list of such anomalously mapped paired-end reads, SVDetect uses a sliding-window strategy to identify all groups of pairs sharing a similar genomic location. The reference genome is divided into overlapping windows of fixed size, and each pair of windows can possibly form a link if at least one pair anchors them by its ends.

## Protocol

The first step in SVDetect is to regroup all pairs that are suspected to originate from the same SV.
The input consists of paired-ends mapped to the reference genome, and the output will contain pairs where either the orientation of pairs is incorrect and/or the distance between them is out of the typical range.

```perl /bioinfo/guests/nriddifo/bin/BAM_preprocessingPairs.pl <sample.sorted.bam>```

and for the reference: 

```perl /bioinfo/guests/nriddifo/bin/BAM_preprocessingPairs.pl <reference.sorted.bam>```

This will create a file containing the aberrant reads for each .bam file: `sample.sorted.ab.bam` and `reference.sorted.bam`, which we then point to in config files.

We also need to create a file called `genome.len` with the number, name and length of each chromosome. This information can be found in the first two columns of the genome.fa.fai file (created by `samtools faidx genome.fa`). `genome.len` must have the following format: 


```
[Chr #]\t[Chr name]\t[Length]
1	2L	3000000
2	2R	4000000
3	3L	5000000
...
```

We need to make a config file for both the sample (tumour) and reference samples that will be used for each step in the analysis. The [SVDetect manual](http://svdetect.sourceforge.net/Site/Manual.html) contains a thorough description of the options for each block 

sample config example:
 

``` bash
<general>
input_format=bam 
sv_type=all
mates_orientation=RF
# Change the next 5 lines accordingly
read1_length=125
read2_length=125
mates_file=/path/to/sample.sorted.ab.bam										
cmap_file=/path/to/genome.len													
output_dir=/path/to/results														
tmp_dir=tmp
num_threads=2
</general>

<detection>
split_mate_file=1
window_size=6541
step_length=1635
</detection>

<filtering>
split_link_file=0
strand_filtering=1
order_filtering=1
insert_size_filtering=1
nb_pairs_threshold=2
nb_pairs_order_threshold=2
indel_sigma_threshold=3
dup_sigma_threshold=2
singleton_sigma_threshold=4
final_score_threshold=0.8
mu_length=2941
sigma_length=233
</filtering>

<bed>
<colorcode>
190,190,190 = 1,2
0,0,0       = 3,3
0,0,255     = 4,4
0,255,0     = 5,5
153,50,205  = 6,7
255,140,0  = 8,10
255,0,0     = 11,10000
</colorcode>
</bed>

<compare>
# Change the next 3 lines accordingly
list_samples=sample.sorted,reference.sorted
list_read_lengths=125-125,125-125
file_suffix=.ab.bam.all.links.filtered
min_overlap=0.05
same_sv_type=1
circos_output=1
bed_output=1
sv_output=1
</compare>
```

A similar file needs to be created for the reference sample. 



















Now, create a bash script to tie them all together: 

```{bash}
#!/bin/bash

WORK_DIR=/path/to/config/files

#Generation and filtering of links from the sample data
#SVDetect linking filtering -conf ${WORK_DIR}/sample.sv.conf

#Generation and filtering of links from the reference data
#SVDetect linking filtering -conf ${WORK_DIR}/reference.sv.conf

#Comparison of links between the two datasets
#SVDetect links2compare -conf ${WORK_DIR}/sample.sv.conf

# Convert to bed format
#SVDetect links2bed -conf ${WORK_DIR}/sample.sv.conf

#Calculation of depth-of-coverage log-ratios
SVDetect cnv ratio2circos ratio2bedgraph -conf ${WORK_DIR}/sample.cnv.conf

#Visualization of filtered links and copy-number profiles in Circos
#circos -conf circos/sample.circos.conf
```





