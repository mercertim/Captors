# Captors
Analysis scripts for Captor pipeline

## Captor analysis pipeline
These scripts a required for the analysis of ONT data using Captors. Demo data (fastq with 2000 entries) and reference files for use with these scripts are contained in this repository. Input and output directories need to be specified by the user along with reference files according to the genome of interest (here we used the Metasequin sequences).

## System requirements and installation guide
The Captor analysis is performed by third party software and does not require installation beyond the below dependencies. RUnning the demo should take <1min. We recommend using a high performance computing (HPC) cluster for the preprocessing step. The preprocessing has been tested on a HPC 64-bit running CentOS Linux release 7.9.2009 with python (3.6.7) and R (4.0.2).

Dependencies command-line:
fastp (v0.20.0) https://github.com/OpenGene/fastp
minimap2 (v2.17-r941) https://github.com/lh3/minimap2
samtools (using htslib v1.9) http://www.htslib.org/
bamtools (v2.2.3) https://github.com/pezmaster31/bamtools

Dependencies Python:
pip install pysamstats (v1.1.2)

Downstream analysis in R can be run on a standard laptop computer and was tested on an Apple Macbook 16G RAM, 500GB memory.

## Demo
A demo dataset for preprocessing containing a Captor library is included.

git clone https://github.com/mercertim/Captors.git

cd Captors

chmod 755 Captor_MetaSequin_perBase_processing.sh

and run the ./Captor_MetaSequin_perBase_processing.sh

The expected output can be found in the file demo_data_captor.bam.bed.tsv

## Instructions for use
NOTE: Generalised scipts for downstream analysis will need to be customised depending on your sample genome. The pipeline is run as follows:
1. Sample preprocessing to determine basewise errors in Captore sequences and target genome using the Captor_MetaSequin_perBase_processing.sh script (BASH and python).

2. Calculate Captor and sample counts using the MetaSequin_ONT_count.R and Captor_ONT_count.R scripts (R).

3. Sample count normalisation using Captors and RUVg the Captor_MetaSequin_normalisation.R script (R).
