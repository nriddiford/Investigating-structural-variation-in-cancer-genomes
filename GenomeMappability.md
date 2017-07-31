# Calculating Genome Mappability with gemtools

## Download gemtools

wget http://barnaserver.com/gemtools/releases/GEMTools-static-i3-1.7.1.tar.gz

# Run gemtools in bash script

```
#PBS -l nodes=1:ppn=16
#PBS -l walltime=1:00:00
#PBS -l mem=2GB
#PBS -j oe /data/kdi_prod/project_result/948/01.00/Analysis/Genomes/GEM/log

genome=/data/kdi_prod/project_result/948/01.00/Analysis/Genomes/Dmel_6/dmel_6.12.fa
out_dir=/data/kdi_prod/project_result/948/01.00/Analysis/Genomes/GEM

export PATH=$HOME/Modules/gemtools-1.7.1-i3/bin:$PATH

log=$out_dir/log
logfile=$log/indexer.log

`> $logfile`

exec &> $log/indexer.runlog.txt

gem-indexer -T 16 -i $genome -o $out_dir/dmel6

gem-mappability -I $out_dir/dmel6.gem -l 200 -o $out_dir/dmel6_mappability_200.gem

gem-2-wig -I $out_dir/dmel6.gem -i $out_dir/dmel6_mappability_200.gem.mappability -o $out_dir/dmel6_mappability

# Convert wig to bed

perl $out_dir/wig2bed.pl $out_dir/dmel6_mappability.wig > $out_dir/dmel6_mappable.bed

# Make 'genome file' from 'mappability.sizes' file
genome_index=/data/kdi_prod/project_result/948/01.00/Analysis/Genomes/Dmel_6/dmel_6.12.fa.fai
awk -v OFS='\t' {'print $1,$2'} $genome_index > $out_dir/dmel6_genome.txt

# make .bed file of regions with low mappability

bedtools complement -i $out_dir/dmel6_mappable.bed -g $out_dir/dmel6_genome.txt > $out_dir/dmel6_unmappable.bed
```
