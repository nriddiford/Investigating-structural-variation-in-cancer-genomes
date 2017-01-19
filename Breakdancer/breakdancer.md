# Tutorial to detect structural variation using breakdancer

# Table of Contents
* [About the tool](#about-the-tool)
* [Input](#input)
* [Run CNV-Seq](#cnv-seq)
* [Output](#output)

* [To do](#to-do)

# About the tool

CNV-Seq allows the detection of copy number variation using NGS data.

# Input

In order to reduce the noise around repeat regions in the genome, first filter reads with a mapping quality > 4: 

`samtools view -b -q 4 sample.bam > sample.qfilt.bam`
