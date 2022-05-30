# Captors
Analysis scripts for Captor pipeline

## Captor analysis pipeline
These scripts a required for the analysis of ONT data using Captors. Data files for use with these scripts are contained in this repository. Input and output directories need to be specified by the user along with reference files according to the genome of interest (here we used the Metasequin sequences).

## System requirements

## Installation guide

## Demo

## Instructions for use
The pipeline is run as follows:
1. Captor preprocessing to determine based wise error rate in the Captor sequences using the Captor_perPore_processing.sh script (BASH and python).

2. Captor preprocessing to determine pore wise error rate in the Captor sequences using the Captor_perPore_processing.sh script (BASH and python).

3. Sample preprocessing to determine basewise errors in target genome using the Captor_MetaSequin_perBase_processing.sh script (BASH and python).

4. Caculate Captor and sample counts using the MetaSequin_ONT_count.R and Captor_ONT_count.R scripts (R).

5. Sample count normalisation using Captors and RUVg the Captor_MetaSequin_normalisation.R script (R).
