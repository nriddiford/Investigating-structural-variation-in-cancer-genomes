# Pipeline

# Table of Contents
* [SVDetector](#svdetector)
  * [About the tool](#about-the-tool)
  * [Protocol](#protocol)

# SVDetector

### About the tool



### Protocol







```fastq_quality_trimmer -t 30 -l 75 -i <input> -o <output>```


```{java}
java -classpath /path/to/Trimmomatic/trimmomatic-0.25.jar org.usadellab.trimmomatic.TrimmomaticPE
-threads 12 \
-phred33 \
<pe_1> <pe_2> <paired_output_1> <unpaired_output_1> <paired_output_2> <unpaired_output_2> \
ILLUMINACLIP:<filter_set> \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:12 MINLEN:36
```


