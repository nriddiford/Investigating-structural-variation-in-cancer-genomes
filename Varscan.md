# Tutorial to detect single nucleotide variation using Varscan2

# Table of Contents
* [About the tool](#about-the-tool)
* [Input](#input)
* [Output](#output)
* [To do](#to-do)

# About the tool


# Input

Varscan takes output from `samtools mplieup`. Easiest to pipe these together: 

`samtools mpileup -f genome.fa sample1.bam -o sample1.mpileup`

This produces a pileup file with per-base information on SNPs:

| chromosome | coordinate | Ref | read_cov | read_bases | base_qualities |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 2L | 228 | T | 4 | .... | BB36 |
| 2L | 229 | G | 4 | .... | BBFF |
| 2L | 230 | A | 5 | ....^!, | BBFFE |
| 2L | 231 | T | 5 | ...., | FFFFF |

Line 1 shows that position 228 on Chr2L is covered by 4 reads, all matching (`.`) the ref nucleotide `T` 


`mpileup2snp sample1.mpileup --min-coverage 25 --min-reads2 4 --min-var-freq 0.1 --p-value 0.05 --output-vcf 1 > sample1.vcf


# Output

Runnig 

| Chrom | Position | Ref | Cons | Reads1 | Reads2 | VarFreq | Strands1 | Strands2 | Qual1 | Qual2 | Pvalue | MapQual1 | MapQual2 | Reads1Plus | Reads1Minus | Reads2Plus | Reads2Minus | VarAllele |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 2L | 5355 | C | Y | 21 | 2 | 8.7% | 2 | 1 | 46 | 69 | 0.98 | 1 | 1 | 14 | 7 | 2 | 0 | T |
| 2L | 5372 | T | W | 6 | 9 | 60% | 2 | 2 | 60 | 44 | 0.98 | 1 | 1 | 5 | 1 | 8 | 1 | A |
| 2L | 5390 | T | W | 13 | 11 | 45.83% | 2 | 2 | 47 | 40 | 0.98 | 1 | 1 | 12 | 1 | 10 | 1 | A |

