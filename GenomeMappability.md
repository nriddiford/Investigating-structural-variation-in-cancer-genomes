# Calculating Genome Mappability with gemtools

## Download gemtools

wget http://barnaserver.com/gemtools/releases/GEMTools-static-i3-1.7.1.tar.gz

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

gem-mappability -I $out_dir/dmel6.gem -l 250 -o $out_dir/mappability_250.gem

gem-2-wig -I dmel6.gem -i mappability_100.gem.mappability -o dmel6_mappability
```
