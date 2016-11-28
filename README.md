# Tools for discovering structural variants in WGS data

# Protocols:
* **SVP callers**
  * [VCFTools](VCFtools.md)
  * [Varscan](Varscan.md)
* **SV callers**
  * [SVDetect](SVDetect.md)
* **CNV callers**
  * [CNV-Seq](CNV-Seq.md)

# Table of contents
* [Background](#background)
* [Overview of detection strategies](#overview-of-detection-strategies)


# Background
Most of the following heavily borrows from [Tattini et al (2015)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4479793/)

Structural variants (SVs) are genomic rearrangements generally affecting more then 50 bp. These can include balanced (where a portion of the genome is moved with no net change in bps) and unbalanced (resulting in an increase/decrease in bps) events  SVs include:
* Deletions
* Insertions
* Inversions
* Mobile-element transpositions
* Translocations
* Tandem repeats
* Copy number variants (CNVs)

# Overview of detection strategies

So far, the NGS-based CNV detection methods can be categorized into five different strategies, including:
* Paired-end mapping (PEM)
* Split read (SR)
* Read depth (RD)
* De novo assembly of a genome (AS)
* Combination of the above approaches (CB)

