# circos plots

Run `run_circos.sh` with the follwing files in working directory:

 o CNV-Seq ['raw_out']: 'HUM-1.tagged.SC.hits.filt-vs-HUM-3.tagged.SC.hits.filt.window-20000.minw-4.cnv'
 o CNV-Seq ['calls']: 'HUM-1_cnvs.txt'
 o Freec ['calls']:'HUM-1.TEx.q4.sort.nodups.RG.bam_CNVs'
 o Lumpy ['Lumpy_output']: 'HUM-1.tagged.SC.lumpy.gt.filtered.vcf'
 o Delly ['Delly_output']:'HUM-1.tagged.SC.gt.delly.vcf'
 o G4s: ['G4s']: 'HUM-1.G4s.bed'
 
`run_circos.sh` will check that each file exists, parse the relevant information and feed it to circos for plotting. 

