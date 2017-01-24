# Protocol to detect copy number variation using Control-FREEC

# Table of Contents
* [About the tool](#about-the-tool)
* [Input](#input)
* [Run Control-Freec](#run-control-freec)
* [Output](#output)
* [Visualising CNVs in IGV](#visualising-cnvs-in-igv)

# About the tool
Control-FREEC is a tool for detection of copy-number changes and allelic imbalances (including LOH) 

# Input
Control-Freec has a very good [manual](http://boevalab.com/FREEC/tutorial.html#CONFIG) that thoroughly explains options that can be set in config file.

Each sample needs an individual config file (start by editing `config_WGS.txt`).

For example see my [example.config.txt](files/example.config.txt) for config for WGS mapped to Drosophila 6.12

# Run Control-Freec

The main script is then called with the config file:

`$ src/./freec -conf example.config.txt`

# Output

Files are written to directory specified in config file. See [manual](http://boevalab.com/FREEC/tutorial.html#OUTPUT) for description of output files.


# Visualising CNVs in IGV

To create an IGV compatible track for each sample, run [cf2gff.pl](script/cf2gff.pl):
