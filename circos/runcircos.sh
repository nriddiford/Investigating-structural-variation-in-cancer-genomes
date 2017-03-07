#!/usr/bin/bash



function usage() {
    echo "
Program: run.circos.sh

This script will look for the following files in the current direrctory, and exit if any one is missing. E.g., for HUM-1. ALL the following files must exist:
	
 o CNV-Seq ['raw_out']: 'HUM-1.tagged.SC.hits.filt-vs-HUM-3.tagged.SC.hits.filt.window-20000.minw-4.cnv'
 o CNV-Seq ['calls']: 'HUM-1_cnvs.txt'
 o Freec ['calls']:'HUM-1.TEx.q4.sort.nodups.RG.bam_CNVs'
 o Lumpy ['Lumpy_output']: 'HUM-1.tagged.SC.lumpy.gt.filtered.vcf'
 o Delly ['Delly_output']:'HUM-1.tagged.SC.gt.delly.vcf'
 o G4s: ['G4s']: 'HUM-1.G4s.bed'
"
}


echo """
################
### CNV-Seq ####
################
"""
cnv_seq_raw_out=$(ls *.cnv | head -n1)

id=$(echo $cnv_seq_raw_out | cut -d '.' -f1 )


cnv_seq_calls=$(ls *_cnvs.txt)

if [[ -f "$cnv_seq_raw_out" && $cnv_seq_calls ]];
then
	echo "Processing CNV-Seq data for sample: '$id'..."

	echo "Reading CNV-Seq raw counts file: '$cnv_seq_raw_out'"
	
	echo "Reading CNV-Seq calls file: '$cnv_seq_calls'"
	
	perl /Users/Nick_curie/Desktop/circos/script/cnv-seq2circos.pl $cnv_seq_raw_out $cnv_seq_calls
else
	echo "CNV-Seq file missing. Aborting"
	usage
	exit 1
fi


echo

echo """
##############
### Freec ####
##############
"""

echo "Processing Freec data"

freec_calls=$(ls *_CNVs | head -n1)

if [[ -f "$freec_calls" ]];
then
	echo "Reading Freec calls file: '$freec_calls'"
	
	perl /Users/Nick_curie/Desktop/circos/script/freec2circos.pl $freec_calls
else
	echo "Freec file '_CNVs' file missing. Aborting"
	usage
	exit 1
fi




echo

echo """
##############
### Lumpy ####
##############
"""

lumpy_filter_set='lumpy_filter_set.txt'

lumpy_characterised_svs='lumpy_characterised_svs.txt'

echo "Processing Lumpy data"

if [ ! -f "$lumpy_filter_set" ]
then
	echo "No filter set provided"
	echo "Creating file '$lumpy_filter_set'"
	touch $lumpy_filter_set
else
	echo "Filtering lumpy calls with '$lumpy_filter_set'"
fi

if [ ! -f "$lumpy_characterised_svs" ]
then
	echo "No characterised variants file provided"
	echo "Creating file '$lumpy_characterised_svs'"
	touch $lumpy_characterised_svs
else 
	echo "Annotating characterised variants with '$lumpy_characterised_svs'"	
fi

lumpy_output=$(ls *.vcf | grep "lumpy")

if [[ -f "$lumpy_output" ]];
then
	echo "Reading Lumpy structural variant calls: '$lumpy_output'"
	
	perl /Users/Nick_curie/Desktop/circos/script/lumpy2circos.pl -v $lumpy_output -d $lumpy_characterised_svs -f $lumpy_filter_set
else
	echo "Lumpy output missing. Aborting"
	usage
	exit 1
fi


echo 

echo """
##############
### Delly ####
##############
"""

delly_filter_set='delly_filter_set.txt'

delly_characterised_svs='delly_characterised_svs.txt'

echo "Processing Delly data"

if [ ! -f "$delly_filter_set" ]
then
	echo "No filter set provided"
	echo "Creating file '$delly_filter_set'"
	touch $delly_filter_set
else
	echo "Filtering lumpy calls with '$delly_filter_set'"
fi

if [ ! -f "$delly_characterised_svs" ]
then
	echo "No characterised variants file provided"
	echo "Creating file '$delly_characterised_svs'"
	touch $delly_characterised_svs
else 
	echo "Annotating characterised variants with '$delly_characterised_svs'"	
fi


delly_output=$(ls *.vcf | grep "delly")

if [[ -f "$delly_output" ]];
then
	echo "Reading Delly structural variant calls: '$delly_output'"
	
	perl /Users/Nick_curie/Desktop/circos/script/delly2circos.pl -v $delly_output -d $delly_characterised_svs -f $delly_filter_set
else
	echo "Delly output missing. Aborting"
	usage
	exit 1
fi

echo

echo """
###################
### Heatmapper ####
###################
"""

echo "Processed all files for $id"

echo "Making heatmap for all variants in 1mb bins"

perl /Users/Nick_curie/Desktop/circos/script/SV_heatmapper.pl



g4s=$(ls *G4s.bed)

if [[ -f "$g4s" ]];
then
	echo "Making heatmap for all G4s in 1mb bins"
	perl /Users/Nick_curie/Desktop/circos/script/SV_heatmapper.pl -f $g4s
else
	echo "G4s bed file missing. Aborting"
	usage
	exit 1
fi


source activate circos


echo """
##############
### Circos ###
##############
"""
echo 

circos && open circos.svg



