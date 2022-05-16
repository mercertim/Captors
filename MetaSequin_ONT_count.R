#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Create tables with Metasequin count per sample
# --------------------------------------------------------------------------

# Set working directory
setwd("path/to/output")

# Load libraries
library('stringr')
library('dplyr')
library('tidyr')
library('ggplot2')
library('rtracklayer')
library('Gviz')
library('ggsci')
library('runner')
library('matrixStats')
library('spgs')
library('ggformula')
library('seqinr')
library('plyr')
library('reshape2')
library('GenomicAlignments')
library('sarlacc')
library('Biostrings')

# Set colour palette
npg_cols <- pal_npg("nrc")(7)

# MetaSequin concentration
sequin_conc <- read.table("MetaSequin_concentration.tsv", header = TRUE)

# List datasets to loop script over
data_sets <- c("Sample_1", "Sample_2", "Sample_3")

# Store kmer analysis for each dataset in lists
sampleList <- vector("list", length(data_sets))

for (samp in data_sets) {
  # Calculate error rates for each base (SUBS and INDELS) for each vector individually
  # --------------------------------------------------------------------------
  
  # Load sequencing error rates
  ONTDNA_pileup_total <- read.table(paste0(samp, ".bam.bed.tsv"), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote="")
  
  # Find reference mapping rates
  get_correct <- function(df){
    bases <- c('A', 'T', 'G', 'C')
    ref <- df['REF_NT']
    correct <- bases[bases == ref]
    sum(as.numeric(df[correct]))
  }
  
  # Number of reasample spanning base
  ONTDNA_pileup_total$Depth <- ONTDNA_pileup_total$Ins + ONTDNA_pileup_total$Del + ONTDNA_pileup_total$Coverage
  
  # Number of reference reasample per base
  ONTDNA_pileup_total$REF_counts <- apply(ONTDNA_pileup_total, 1, get_correct)
  
  # Number rate of reference mapping as a fraction of total spanning reasample
  ONTDNA_pileup_total$REF_rate <- ONTDNA_pileup_total$REF_counts / ONTDNA_pileup_total$Depth
  
  # Find substitution rates
  get_subs <- function(df){
    bases <- c('A', 'T', 'G', 'C')
    ref <- df['REF_NT']
    subs <- bases[!(bases == ref)]
    sum(as.numeric(df[subs]))
  }
  
  # Number and rate of substitutions per base
  ONTDNA_pileup_total$SUB_counts <- apply(ONTDNA_pileup_total, 1, get_subs)
  ONTDNA_pileup_total$SUB_rate <- ONTDNA_pileup_total$SUB_counts / ONTDNA_pileup_total$Depth
  
  # Number and rate of INDELs per base
  ONTDNA_pileup_total$Ins_rate <- ONTDNA_pileup_total$Ins / ONTDNA_pileup_total$Depth
  ONTDNA_pileup_total$Del_rate <- ONTDNA_pileup_total$Del / ONTDNA_pileup_total$Depth
  ONTDNA_pileup_total$INDEL_counts <- ONTDNA_pileup_total$Ins + ONTDNA_pileup_total$Del
  ONTDNA_pileup_total$INDEL_rate <- ONTDNA_pileup_total$INDEL_counts / ONTDNA_pileup_total$Depth
  
  # Number and rate of errors per base
  ONTDNA_pileup_total$Error_counts <- ONTDNA_pileup_total$SUB_counts + ONTDNA_pileup_total$INDEL_counts
  ONTDNA_pileup_total$Error_rate <- ONTDNA_pileup_total$Error_counts / ONTDNA_pileup_total$Depth
  
  # Calculate GC content indepenedent and Homopolymer content of 30mers
  # --------------------------------------------------------------------------
  
  for(sequin in unique(ONTDNA_pileup_total$Chrom)) {
    
    # Subset to each genome
    ONTDNA_pileup <- ONTDNA_pileup_total[ONTDNA_pileup_total$Chrom == sequin, ]
    
    # Calculate all unique 30mers in Metasequin sequences
    adaptors_seq_F <- ONTDNA_pileup$REF_NT
    
    KD_k <- function(x){paste(x, collapse = "")}
    adaptors_seq_F_30mers <- KD_k(adaptors_seq_F)
    
    # Make ranges of 30mers in Metasequin sequences
    k_total_granges <- GRanges(seqnames = sequin, ranges = IRanges(start = 1:length(adaptors_seq_F_30mers), width = length(adaptors_seq_F_30mers)), strand = "+", SEQUENCE = adaptors_seq_F_30mers)
    
    # Calculate error rates
    k_total_granges$Error_freq <- sum(extractList(ONTDNA_pileup$Error_counts, k_total_granges@ranges))
    k_total_granges$Error_mean <- mean(extractList(ONTDNA_pileup$Error_rate, k_total_granges@ranges))
    k_total_granges$Error_sd <- sd(extractList(ONTDNA_pileup$Error_rate, k_total_granges@ranges))
    k_total_granges$Error_upper <- k_total_granges$Error_mean + k_total_granges$Error_sd
    k_total_granges$Error_lower <- k_total_granges$Error_mean - k_total_granges$Error_sd
    
    # Calulcate coverage (non-indel spanning reads)
    k_total_granges$Coverage_total <- sum(extractList(ONTDNA_pileup$Coverage, k_total_granges@ranges))
    k_total_granges$Coverage_mean <- mean(extractList(ONTDNA_pileup$Coverage, k_total_granges@ranges))
    k_total_granges$Coverage_max <- max(extractList(ONTDNA_pileup$Coverage, k_total_granges@ranges))
    
    # Calculate depth (all spanning reads)
    k_total_granges$Depth_total <- sum(extractList(ONTDNA_pileup$Depth, k_total_granges@ranges))
    k_total_granges$Depth_mean <- mean(extractList(ONTDNA_pileup$Depth, k_total_granges@ranges))
    k_total_granges$Depth_max <- max(extractList(ONTDNA_pileup$Depth, k_total_granges@ranges))
    
    
    # Add to list for each MetaSequin per sample
    sampleList[[paste0(samp, "_", sequin)]] <- k_total_granges
  }
}

# Make info tables
# --------------------------------------------------------------------------

# Make dataframe of sequin coverage
for(i in names(sampleList)){
  tmp <- sampleList[[i]]
  sampleName <- i
  name <- as.character(tmp@seqnames@values)
  Depth <- tmp$Depth_mean
  Depth_max <- tmp$Depth_max
  Errors <- tmp$Error_mean
  
  if(i == names(sampleList)[1]){
    metaseqin_depth <- data.frame("samp" = sampleName, 'Adaptor' = name, 'Depth' = Depth, 'Depth_max' = Depth_max, 'Error' = Errors)
  } else {
    metaseqin_depth <- rbind(metaseqin_depth, data.frame("samp" = sampleName, 'Adaptor' = name, 'Depth' = Depth, 'Depth_max' = Depth_max, 'Error' = Errors))
  }
}

# Add metasequin concentration based on the mixture used in sample
idx <- match(metaseqin_depth$Adaptor, sequin_conc$ID)
metaseqin_depth$MIX_A <- sequin_conc$MIX_A [idx]
metaseqin_depth$MIX_B <- sequin_conc$MIX_B [idx]
metaseqin_depth <- metaseqin_depth[!(metaseqin_depth$Adaptor == 'SynX' | metaseqin_depth$Adaptor == 'PhiX'), ]
metaseqin_depth <- metaseqin_depth[grepl("MG|ML", metaseqin_depth$Adaptor), ]
metaseqin_depth$Mix <- ifelse(grepl("A[1-3]", metaseqin_depth$samp), "A", "B")

for(i in 1:length(metaseqin_depth$Mix)){
  if(i == 1){
    conc <- c()
  }
  
  if(metaseqin_depth$Mix[i] == "A"){
    conc <- c(conc, metaseqin_depth$MIX_A[i])
  } else {
    conc <- c(conc, metaseqin_depth$MIX_B[i])
  }
}
metaseqin_depth$CONC <- conc

#  Write tables with Metasequin depth (mean and max)
metaseqin_depth_mean <- dcast(metaseqin_depth, Adaptor ~ samp, value.var = "Depth")
idx <- match(metaseqin_depth_mean$Adaptor, sequin_conc$ID)
metaseqin_depth_mean$MIX_A <- sequin_conc$MIX_A [idx]
metaseqin_depth_mean$MIX_B <- sequin_conc$MIX_B [idx]
metaseqin_depth_mean[is.na(metaseqin_depth_mean)] <- 0
write.csv(metaseqin_depth_mean, "Metasequin_coverge_by_sample.csv", row.names = FALSE)



