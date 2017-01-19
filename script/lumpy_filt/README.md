# lumpy_filt.pl

This program filters the output from LUMPY.

# LUMPY output



 when run on a tumour:normal pair as follows:

o Map reads with BWA (with RG info)

o Run LUMPY express:
	 lumpyexpress \
	    -B tumour.bam,normal.bam \
	    -S tumour.split.sort.bam,normal.split.sort.bam \
    	    -D tumour.discordant.sort.bam,normal.discordant.sort.bam \
    	    -o tumour.normal.vcf

o Genotype output with SVTyper:
	 svtools genotype \
	    -i tumour.normal.vcf
	    -o tumour.normal.gt.vcf
	    -B tumour.bam,normal.bam

Run as perl lumpy_filt.pl -h to see usage statement

Nick Riddiford 2017