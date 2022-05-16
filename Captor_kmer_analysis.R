#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Analysis of Captor 6-mers
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
library('ggpubr')
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

# List datasets to loop script over
data_sets <- c("Sample_1", "Sample_2", "Sample_3")

# Read depth information
ONT_adaptors <- read.csv("Captors_concentration.csv", header = TRUE)

# Store 6-mer analysis for each dataset in lists
sampleList <- vector("list", length(data_sets))

for (samp in data_sets) {
  # Calculate error rates for each base (SUBS and INDELS) for each vector individually
  # --------------------------------------------------------------------------
  
  # Load sequencing error rates from perBase processing script
  ONTDNA_pileup_total <- read.table(paste0(samp, "_adaptor.bam.bed.tsv"), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote="")
  
  # Subset to Captor variable region
  ONTDNA_pileup_total <- ONTDNA_pileup_total[ONTDNA_pileup_total$Number >= 31 & ONTDNA_pileup_total$Number <= 60,]
  
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
  
  # Find substitution rates per base
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
  
  # Calculate GC content and Homopolymer content of 6-mers
  # --------------------------------------------------------------------------
  
  # Loop for each Captor
  for(captor in unique(ONTDNA_pileup_total$Chrom)) {
    
    # Subset to each Captor
    ONTDNA_pileup <- ONTDNA_pileup_total[ONTDNA_pileup_total$Chrom == captor, ]
    
    # Calculate all unique 6-mer in Captor
    Captor_seq_F <- ONTDNA_pileup$REF_NT
    
    KD_k <- function(x){paste(x, collapse = "")}
    Captor_seq_F_6mers <- runner(x = Captor_seq_F, k = 6, f = KD_k)
    Captor_seq_total_6mers <- Captor_seq_F_6mers[str_length(Captor_seq_F_6mers) == 6]
    
    # Make ranges of 6-mers in Captors sequence
    k_total_granges <- GRanges(seqnames = captor, ranges = IRanges(start = 1:length(Captor_seq_total_6mers), width = 6), strand = "+", SEQUENCE = Captor_seq_total_6mers)
    
    # Calculate GC content
    GC_k <- function(x){(str_count(x, pattern = 'C') + str_count(x, pattern = 'G')) / 6 *100}
    k_total_granges$GC_pct <- GC_k(k_total_granges$SEQUENCE)
    
    # Calculate homopolymers in 6-mers
    hp_max <- data.frame(k_total_granges$SEQUENCE)
    hp_max$HP2 <- grepl("AA|TT|GG|CC", k_total_granges$SEQUENCE)
    hp_max$HP3 <- grepl("AAA|TTT|GGG|CCC", k_total_granges$SEQUENCE)
    hp_max$HP4 <- grepl("AAAA|TTTT|GGGG|CCCC", k_total_granges$SEQUENCE)
    hp_max$HP5 <- grepl("AAAAA|TTTTT|GGGGG|CCCCC", k_total_granges$SEQUENCE)
    hp_max$HP6 <- grepl("AAAAAA|TTTTTT|GGGGGG|CCCCCC", k_total_granges$SEQUENCE)
    k_total_granges$HP_length <- rowSums(hp_max[,2:6]) + 1
    
    # Add Captors concentration
    idx <- match(k_total_granges@seqnames@values, ONT_adaptors$Captor_ID)
    k_total_granges$Captor_conc <- ONT_adaptors$CONC [idx]
    
    # Calculate error rates
    k_total_granges$Error_freq <- sum(extractList(ONTDNA_pileup$Error_counts, k_total_granges@ranges))
    k_total_granges$Error_mean <- mean(extractList(ONTDNA_pileup$Error_rate, k_total_granges@ranges))
    k_total_granges$Error_sd <- sd(extractList(ONTDNA_pileup$Error_rate, k_total_granges@ranges))
    k_total_granges$Error_upper <- k_total_granges$Error_mean + k_total_granges$Error_sd
    k_total_granges$Error_lower <- k_total_granges$Error_mean - k_total_granges$Error_sd
    
    k_total_granges$SUB_freq <- sum(extractList(ONTDNA_pileup$SUB_counts, k_total_granges@ranges))
    k_total_granges$SUB_mean <- mean(extractList(ONTDNA_pileup$SUB_rate, k_total_granges@ranges))
    k_total_granges$SUB_sd <- sd(extractList(ONTDNA_pileup$SUB_rate, k_total_granges@ranges))
    k_total_granges$SUB_upper <- k_total_granges$SUB_mean + k_total_granges$SUB_sd
    k_total_granges$SUB_lower <- k_total_granges$SUB_mean - k_total_granges$SUB_sd
    
    k_total_granges$INDEL_freq <- sum(extractList(ONTDNA_pileup$INDEL_counts, k_total_granges@ranges))
    k_total_granges$INDEL_mean <- mean(extractList(ONTDNA_pileup$INDEL_rate, k_total_granges@ranges))
    k_total_granges$INDEL_sd <- sd(extractList(ONTDNA_pileup$INDEL_rate, k_total_granges@ranges))
    k_total_granges$INDEL_upper <- k_total_granges$INDEL_mean + k_total_granges$INDEL_sd
    k_total_granges$INDEL_lower <- k_total_granges$INDEL_mean - k_total_granges$INDEL_sd
    
    k_total_granges$INS_freq <- sum(extractList(ONTDNA_pileup$Ins, k_total_granges@ranges))
    k_total_granges$INS_mean <- mean(extractList(ONTDNA_pileup$Ins_rate, k_total_granges@ranges))
    k_total_granges$INS_sd <- sd(extractList(ONTDNA_pileup$Ins_rate, k_total_granges@ranges))
    k_total_granges$INS_upper <- k_total_granges$INS_mean + k_total_granges$INS_sd
    k_total_granges$INS_lower <- k_total_granges$INS_mean - k_total_granges$INS_sd
    
    k_total_granges$DEL_freq <- sum(extractList(ONTDNA_pileup$Del, k_total_granges@ranges))
    k_total_granges$DEL_mean <- mean(extractList(ONTDNA_pileup$Del_rate, k_total_granges@ranges))
    k_total_granges$DEL_sd <- sd(extractList(ONTDNA_pileup$Del_rate, k_total_granges@ranges))
    k_total_granges$DEL_upper <- k_total_granges$DEL_mean + k_total_granges$DEL_sd
    k_total_granges$DEL_lower <- k_total_granges$DEL_mean - k_total_granges$DEL_sd
    
    k_total_granges$Coverage_total <- sum(extractList(ONTDNA_pileup$Coverage, k_total_granges@ranges))
    k_total_granges$Coverage_mean <- mean(extractList(ONTDNA_pileup$Coverage, k_total_granges@ranges))
    
    k_total_granges$Depth_total <- sum(extractList(ONTDNA_pileup$Depth, k_total_granges@ranges))
    k_total_granges$Depth_mean <- mean(extractList(ONTDNA_pileup$Depth, k_total_granges@ranges))

    # Calculate error rates across 6-mers
    # --------------------------------------------------------------------------
    
    # Calculate mean error rate per 6-mer in each Captors (for when a 6-mer occurs more than once in a Captor)
    kmer_error_freq <- data.frame(k_total_granges@elementMetadata) %>%
      group_by(SEQUENCE) %>%
      dplyr::summarize(Error_rate_mean = mean(Error_mean, na.rm = TRUE),
                       SUB_rate_mean = mean(SUB_mean, na.rm = TRUE),
                       INS_rate_mean = mean(INS_mean, na.rm = TRUE),
                       DEL_rate_mean = mean(DEL_mean, na.rm = TRUE),
                       GC_pct = mean(GC_pct, na.rm = TRUE),
                       HP_length = mean(HP_length, na.rm = TRUE),
                       max_conc = max(Captor_conc))
    
    # Add to list for each Captor per sample
    sampleList[[paste0(samp, "_", captor)]] <- kmer_error_freq
  }
}

# Make a dataframe of 6-mers for each Captor per sample
for(i in grep('samp', names(sampleList), value = TRUE)){
  tmp <- data.frame(sampleList[[i]])
  tmp$Name <- i
  
  if(i == grep('samp', names(sampleList), value = TRUE)[1]){
    Captor_kmers <- tmp
  } else {
    Captor_kmers <- rbind(Captor_kmers, tmp)
  }
}

# Collate errors per 6-mer
Captor_kmer_error_mean <- dcast(data.frame(Captor_kmers), SEQUENCE ~ Name, value.var = "Error_rate_mean", fun.aggregate = mean)

# Add maximum Captor concentration for each 6-mer
max_kmer_conc <- data.frame(Captor_kmers) %>%
  group_by(SEQUENCE) %>%
  dplyr::summarize(max_conc = max(max_conc)) %>%
  select(SEQUENCE, max_conc)

idx <- match(Captor_kmer_error_mean$SEQUENCE, max_kmer_conc$SEQUENCE)
Captor_kmer_error_mean$Max_Captor_concentration <- max_kmer_conc$max_conc [idx]

write.csv(Captor_kmer_error_mean, "6mers_mean_error_rate_per_sample.csv", row.names = FALSE)
