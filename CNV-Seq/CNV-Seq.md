# Tutorial to detect copy number variation using CNV-Seq

# Table of Contents
* [About the tool](#about-the-tool)
* [Input](#input)
  * [Hits files](#hits-file)
* [Run CNV-Seq](#run-cnv-seq)
* [Output](#output)
* [Plotting](#plotting)

# About the tool

CNV-Seq allows the detection of copy number variation using NGS data.

# Input

In order to reduce the noise around repeat regions in the genome, first filter reads with a mapping quality > 4: 

`samtools view -b -q 4 sample.bam > sample.qfilt.bam`

## Hits file
We only need to provide two "hits" files for CNV-Seq to work, one for the sample, and one for the reference. 

These can be extracted as follows: 

`samtools view -F 4 sample.qfilt.bam | perl -lane 'print "$F[2]\t$F[3]"' > sample.qfilt.hits` 

This creates a two column, tab delimited file, with the second column giving the corresponding 1-based leftmost mapping position of a read.

```
2L	1
2L	1
2L	4
2L	4
2L	4
2L	4
2L	4
2L	4
2L	4 
```

Next, to select only fully assembled chromosomes run [filt.pl](script/filt.pl). This will output filtered files eg `sample.qfilt.hits.filt`:


# Run CNV-Seq
Now run the main perl script: 

`cnv-seq.pl --ref ref.qfilt.hits --test sample.qfilt.hits --genome-size 23542271 window-size 10000`

# Output

This produces two files `sample-vs-reference.cnv` and `sample-vs-reference.count`

`sample-vs-reference.count` shows the raw count data for each CNV. 

| chromosome | start | end | test | ref |
|:---:|:---:|:---:|:---:|:---:|
| X | 1 | 363 | 70 | 124 |
| X | 183 | 545 | 82 | 123 |
| X | 365 | 727 | 90 | 115 |

`sample-vs-reference.cnv` contains the stats. 

| chromosome | start | end | test | ref | position | log2 | p.value | cnv | cnv.size | cnv.log2 | cnv.p.value |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| X | 1 | 363 | 70 | 124 | 182 | -0.479319752689881 | 0.0462467378993667 | 0 | NA | NA | NA |
| X | 183 | 545 | 82 | 123 | 364 | -0.239368959969129 | 0.194296291525077 | 0 | NA | NA | NA |
| X | 365 | 727 | 90 | 115 | 546 | -0.00804341386267303 | 0.488268352520283 | 0 | NA | NA | NA |


# Plotting

To plot and output CNV info for samples, run [cnv_seq_process.sh](script/cnv_seq_process.sh) on `.cnv` files:

`$ bash cnv_seq_process.sh *.cnv`

This uses a modified version of the main script provided with [CNV-Seq](https://github.com/hliang/cnv-seq/blob/master/cnv/R/cnv.R).

The main tweaks:
  * Change the plotting colours/densities
  * Show genes of interest in a closeup plot of notch region [e.g.](files/HUM-7_notch.pdf)
  * Save CNV details for each file parsed for later processing

This will by default produce two plots (one for [notch region](files/HUM-7_notch.pdf), and one for [chromosome X](files/HUM-7_X.pdf)) and a [cnvs.txt](HUM-4_cnvs.txt) file. 
To plot for different chromosomes alter the `chrom` var in the script [cnv_seq_process.sh](script/cnv_seq_process.sh).




