# Captors
Analysis scripts for Captor pipeline

## Captor analysis pipeline
These scripts a required for the analysis of ONT data using Captors. Demo data (fastq with 2000 entries) and reference files for use with these scripts are contained in this repository. Input and output directories need to be specified by the user along with reference files according to the genome of interest (here we used the Metasequin sequences).

## System requirements
The Captor analysis is performed by third party software and does not require installation beyond the below dependencies. We recommend using a high performance computing (HPC) cluster for the preprocessing steps. The preprocessing has been tested on a sungrid engine HPC 64-bit running python (3.6) and R (4.0.2).

Dependencies command-line:
fastp https://github.com/OpenGene/fastp
minimap2 https://github.com/lh3/minimap2
samtools http://www.htslib.org/
bamtools https://github.com/pezmaster31/bamtools

Dependencies Python:
pip install pysamstats

Downstream analysis can be run on a standard laptop computer and was tested on an Apple Macbook 16G RAM, 500GB memory. Dependencies R: 

## Installation guide

## Demo

## Instructions for use
The pipeline is run as follows:
1. Captor preprocessing to determine based wise error rate in the Captor sequences using the Captor_perPore_processing.sh script (BASH and python).

2. Captor preprocessing to determine pore wise error rate in the Captor sequences using the Captor_perPore_processing.sh script (BASH and python).

3. Sample preprocessing to determine basewise errors in target genome using the Captor_MetaSequin_perBase_processing.sh script (BASH and python).

4. Caculate Captor and sample counts using the MetaSequin_ONT_count.R and Captor_ONT_count.R scripts (R).

5. Sample count normalisation using Captors and RUVg the Captor_MetaSequin_normalisation.R script (R).
