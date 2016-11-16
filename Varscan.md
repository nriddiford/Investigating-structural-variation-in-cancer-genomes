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

`varscan mpileup2snp sample1.mpileup --min-coverage 25 --min-reads2 4 --min-var-freq 0.1 --p-value 0.05 --output-vcf 1 > sample1.vcf`

These are quite consservative settings, for a position to be called as an SNP in sample1, it needs to have:
* ≥25× read coverage at that site
* Variant must be supported by at least 4 reads
* Must have an an allele frequency of ≥0.1
* P-value ≤0.05

## somatic

Running varscan with somatic option allows the comparison between two samples: 

`varscan somatic sample1.mpileup control.mpileup S1vsC --min-coverage 20 --tumor-purity 0.75`

This gives two files, one called `S1vsC.snp` and one `S1vsC.indel` with the following collumns: 

1. **Chrom**		chromosome name
2.	**Position**	position (1-based)
3.	**Ref**		reference allele at this position
4.	**Var**		variant allele at this position
5.	**Normal_Reads1**	reads supporting reference allele
6.	**Normal_Reads2**	reads supporting variant allele
7.	**Normal_VarFreq**	frequency of variant allele by read count
8.	**Normal_Gt**	genotype call for Normal sample
9.	**Tumor_Reads1**	reads supporting reference allele
10.	**Tumor_Reads2**	reads supporting variant allele
11.	**Tumor_VarFreq**	frequency of variant allele by read count
12.	**Tumor_Gt**	genotype call for Tumor sample
13.	**Somatic_Status**	status of variant (Germline, Somatic, or LOH)	
14.	**Pvalue**		Significance of variant read count vs. expected baseline error
15.	**Somatic_Pvalue**	Significance of tumor read count vs. normal read count

| chrom | position | ref | var | normal_reads1 | normal_reads2 | normal_var_freq | normal_gt | tumor_reads1 | tumor_reads2 | tumor_var_freq | tumor_gt | somatic_status | variant_p_value | somatic_p_value | tumor_reads1_plus | tumor_reads1_minus | tumor_reads2_plus | tumor_reads2_minus | normal_reads1_plus | normal_reads1_minus | normal_reads2_plus | normal_reads2_minus |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 2L | 5390 | T | A | 26 | 16 | 38.1% | W | 13 | 11 | 45.83% | W | Germline | 2.580689727008415E-10 | 0.36014277004737094 | 12 | 1 | 10 | 1 | 14 | 12 | 9 | 7 |
| 2L | 5403 | C | G | 27 | 13 | 32.5% | S | 15 | 11 | 42.31% | S | Germline | 4.589401497448938E-9 | 0.29105012328732616 | 14 | 1 | 9 | 2 | 13 | 14 | 8 | 5 |
| 2L | 5465 | C | A | 26 | 15 | 36.59% | M | 18 | 9 | 33.33% | M | Germline | 5.034423021681378E-9 | 0.7016421521628455 | 12 | 6 | 9 | 0 | 13 | 13 | 12 | 3 |

And one called `S1vsC.indel`:

| chrom | position | ref | var | normal_reads1 | normal_reads2 | normal_var_freq | normal_gt | tumor_reads1 | tumor_reads2 | tumor_var_freq | tumor_gt | somatic_status | variant_p_value | somatic_p_value | tumor_reads1_plus | tumor_reads1_minus | tumor_reads2_plus | tumor_reads2_minus | normal_reads1_plus | normal_reads1_minus | normal_reads2_plus | normal_reads2_minus |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 2L | 9223 | A | -TGT | 19 | 9 | 32.14% | */-TGT | 8 | 14 | 63.64% | */-TGT | Germline | 1.0 | 0.026264528753106215 | 1 | 7 | 8 | 6 | 11 | 8 | 4 | 5 |
| 2L | 16479 | A | +T | 43 | 13 | 23.21% | */+T | 32 | 11 | 25.58% | */+T | Germline | 1.219798980988551E-8 | 0.4835645873043649 | 21 | 11 | 10 | 1 | 26 | 17 | 13 | 0 |
| 2L | 16721 | C | -A | 20 | 12 | 37.5% | */-A | 19 | 6 | 24% | */-A | Germline | 7.746173604018616E-7 | 0.916488594650541 | 11 | 8 | 3 | 3 | 14 | 6 | 3 | 9 |

# To do

