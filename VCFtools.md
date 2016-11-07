# Calling SNPs with VCFtools

# Table of Contents
* [About the tool](#about-the-tool)
* [Protocol](#protocol)

## About the tool
[VCFtools](https://vcftools.github.io) is a suite of tools that work with .vcf files

## Protocol

VCFtools (and the binary version BCFtools) accept .vcf fils that can generated through [samtools mpileup][http://www.htslib.org/doc/samtools.html] (for example).

.bam files must be in the same directory as their index file

```samtools mpileup -ugf genome.fa sample1_sorted_.bam | bcftools call -vmO z -o sample1_sorted_variants.vcf.gz```

Then the .vcf files need to be indexed: 

```tabix -p vcf sample1_sorted_variants.vcf.gz``` 

We can then generate a .bedfile for viewing in IGV: 

```bedtools genomecov -ibam sample1_sorted_.bam -bg -trackline -trackopts 'sample1_sorted_variants" visibility=1 color=255,30,30' > sample1_sorted_variants.bedgraph```

If we want to find SNPs unique to 

```bcftools isec sample1_sorted_variants.vcf.gz sample2_sorted_variants.vcf.gz sample3_sorted_variants.vcf.gz -p Diff/ -n -1 -c all```

This will create new directory called Diff, containing the results. SNPS unique to sample1 will be in `0000.vcf`. This can then be indexed using IGVtools (load as file and then run IGVtools "toTDF"), and then used as an annotation track.