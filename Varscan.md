# Tutorial to detect single nucleotide variation using Varscan2

# Table of Contents
* [About the tool](#about-the-tool)
* [Input](#input)
  * [mpileup](#mpileup)
* [Run Varscan](#run-varscan)
  * [mpileup2snp](#mpileup2snp)
  * [somatic](#somatic)
* [To do](#to-do)

# About the tool


# Input

Varscan takes output from `samtools mplieup`. E.g:

## Mpileup

`samtools mpileup -f genome.fa sample1.bam -o sample1.mpileup`

This produces a pileup file with per-base information in sample1:

| chromosome | coordinate | Ref | read_cov | read_bases | base_qualities |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 2L | 228 | T | 4 | .... | BB36 |
| 2L | 229 | G | 4 | .... | BBFF |
| 2L | 230 | A | 5 | ....^!, | BBFFE |
| 2L | 231 | T | 5 | ...., | FFFFF |

Row 1 for example shows that the 1-based coordinate 228 on Chr2L is covered by 4 reads, all matching (`.`) the ref nucleotide `T` 

# Run Varscan

## mpileup2snp

Run mpileup2snp to identify SNPs in a particular sample (vs reference genome):

`mpileup2snp sample1.mpileup --min-coverage 25 --min-reads2 4 --min-var-freq 0.1 --p-value 0.05 --output-vcf 1 > sample1.vcf`

These are quite consservative settings, for a position to be called as an SNP in sample1, it needs to have:
* ≥25× read coverage at that site
* Variant must be supported by at least 4 reads
* Must have an an allele frequency of ≥0.1
* P-value ≤0.05

## somatic

Running varscan with somatic option  

Runnig 

| chrom | position | ref | var | normal_reads1 | normal_reads2 | normal_var_freq | normal_gt | tumor_reads1 | tumor_reads2 | tumor_var_freq | tumor_gt | somatic_status | variant_p_value | somatic_p_value | tumor_reads1_plus | tumor_reads1_minus | tumor_reads2_plus | tumor_reads2_minus | normal_reads1_plus | normal_reads1_minus | normal_reads2_plus | normal_reads2_minus |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 2L | 5390 | T | A | 26 | 16 | 38.1% | W | 13 | 11 | 45.83% | W | Germline | 2.580689727008415E-10 | 0.36014277004737094 | 12 | 1 | 10 | 1 | 14 | 12 | 9 | 7 |
| 2L | 5403 | C | G | 27 | 13 | 32.5% | S | 15 | 11 | 42.31% | S | Germline | 4.589401497448938E-9 | 0.29105012328732616 | 14 | 1 | 9 | 2 | 13 | 14 | 8 | 5 |
| 2L | 5465 | C | A | 26 | 15 | 36.59% | M | 18 | 9 | 33.33% | M | Germline | 5.034423021681378E-9 | 0.7016421521628455 | 12 | 6 | 9 | 0 | 13 | 13 | 12 | 3 |

# To do

