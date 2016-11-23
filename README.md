# Tools for discovering structural variants in WGS data

* **SVP callers**
  * [VCFTools](VCFtools.md)
  * [Varscan](Varscan.md)
* **SV callers**
  * [SVDetect](SVDetect.md)
* **CNV callers**
  * [CNV-Seq](CNV-Seq.md)

# Background
Most of the following heavilly borrows from (Tattini et al (2015))[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4479793/]

Structural variants (SVs) are genomic rearrangements generally affecting more then 50 bp. SVs include:
* Deletions
* Insertions
* Inversions
* Mobile-element transpositions
* Translocations
* Tandem repeats
* Copy number variants (CNVs)

Several databases – e.g., the Database of Genomic Variants archive which reports structural variation identified in healthy control samples (DGVa1) – have been created for the collection of SVs data (Lappalainen et al., 2013). Public data resources have been developed with the purpose of supporting the interpretation of clinically relevant variants, e.g., dbVar2, or collecting known disease genes (OMIM3) hit by SVs.

