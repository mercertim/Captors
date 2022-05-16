#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Normalise Metasequin depth based on adaptor counts per sample
# --------------------------------------------------------------------------

# Set adaptive working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/tim_projects/adaptor/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/tim/mer/scott/adaptor/")
  place <- "timmer"
}

# Libraries
library('edgeR')
library('limma')
library('RUVSeq')
library('ggplot2')
library('ggrepel')
library('tidyr')
library('dplyr')
library('ggsci')
library('rtracklayer')
library('GenomicRanges')
library('GenomicFeatures')
library('subSeq')
library('DESeq2')

# Load and prepare data
# --------------------------------------------------------------------------

# Read depth information
sequin_conc <- read.table("/Users/mac/cloudstor/tim_projects/synx/data/reference_files/metasequin_features_corrected.tsv", header = TRUE)

# Load dataframes of counts
adaptors <- read.csv("project_results/adaptor_concentration/Adaptor_coverge_by_sample.csv", header = TRUE, row.names = 1)
metaseqins <- read.csv("project_results/metasequin_concentration/Metasequin_coverge_by_sample.csv", header = TRUE, row.names = 1)
total_counts <- round(rbind(adaptors[,1:6], metaseqins[,1:6]))

# Create sample design matrix
sampleName <- gsub("ONT_SeqMetGen", "Mix", colnames(total_counts))
sampleType <- factor(c(gsub( "[0-9]$", "", sampleName)))
sampleRep <- factor(gsub("[^0-9]", "", sampleName))
sampleTable <- data.frame(sampleName, sampleType, sampleRep)
rownames(sampleTable) <- colnames(total_counts)

# Plot RLE without any normalisation
set <- newSeqExpressionSet(as.matrix(total_counts), phenoData = data.frame(colnames(total_counts), row.names=colnames(total_counts)))
pdf("project_results/metasequin_concentration/RLE_UNnormalised.pdf")
plotRLE(set, outline=FALSE)
dev.off()

# Log-fold change without normalisation
y <- DGEList(counts = as.matrix(total_counts), group=sampleType)
design <- model.matrix(~sampleType)
y <- estimateDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
UNres <- data.frame(topTags(lrt, n = Inf))

# Add expected LFC information
idx <- match(rownames(UNres), sequin_conc$ID)
UNres$MIX_A <- sequin_conc$MIX_A [idx]
UNres$MIX_B <- sequin_conc$MIX_B [idx]
UNres$MIX_A[is.na(UNres$MIX_A)] <- 1
UNres$MIX_B[is.na(UNres$MIX_B)] <- 1
UNres$Expected_LFC <- log2(UNres$MIX_B/UNres$MIX_A)
UNres$Spikes <- grepl("^ONT", rownames(UNres))

# Linear model of observed vs expected
mod <- lm(logFC ~ Expected_LFC, UNres)
eq <- paste0("y = ", format(unname(coef(mod)[1]), digits = 2), " + ", format(unname(coef(mod)[2]), digits = 2), "x")
r2 <- paste0("r2 ", format(summary(mod)$r.squared, digits = 3))
eq_r2 <- paste0(eq, " ", r2)

# Plot observed vs expected
ggplot(UNres, aes(x=Expected_LFC, y=logFC)) +
  annotate(geom = "text", y = 3, x = -3, label = eq_r2) +
  geom_point(aes(color = Spikes), shape = 1) +
  geom_smooth(method='lm', se=TRUE, color = 'blue', linetype = 'dashed') +
  theme_classic() +
  theme(aspect.ratio = 1) +
  xlab("Expected LFC") +
  ylab("Observed LFC") +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed')
ggsave("project_results/metasequin_concentration/MixA_vs_MixB_LFC_UNnormalise_edgeR.pdf")

# Perform TMM normalisation and DE testing
y <- DGEList(counts = as.matrix(total_counts), group=sampleType)
y <- calcNormFactors(y, method="TMM")
design <- model.matrix(~sampleType)
y <- estimateDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
TMMres <- data.frame(topTags(lrt, n = Inf))

# Plot RLE with TMM normalisation
set <- newSeqExpressionSet(as.matrix(total_counts), normalizedCounts = round(cpm(y, normalized.lib.sizes = TRUE)), phenoData = data.frame(colnames(total_counts), row.names=colnames(total_counts)))

pdf("project_results/metasequin_concentration/RLE_TMM_normalised.pdf")
plotRLE(set, outline=FALSE)
dev.off()

# Add expected LFC information
idx <- match(rownames(TMMres), sequin_conc$ID)
TMMres$MIX_A <- sequin_conc$MIX_A [idx]
TMMres$MIX_B <- sequin_conc$MIX_B [idx]
TMMres$MIX_A[is.na(TMMres$MIX_A)] <- 1
TMMres$MIX_B[is.na(TMMres$MIX_B)] <- 1
TMMres$Expected_LFC <- log2(TMMres$MIX_B/TMMres$MIX_A)
TMMres$Spikes <- grepl("^ONT", rownames(TMMres))

# Linear model of observed vs expected
mod <- lm(logFC ~ Expected_LFC, TMMres)
eq <- paste0("y = ", format(unname(coef(mod)[1]), digits = 2), " + ", format(unname(coef(mod)[2]), digits = 2), "x")
r2 <- paste0("r2 ", format(summary(mod)$r.squared, digits = 3))
eq_r2 <- paste0(eq, " ", r2)

# Plot observed vs expected
ggplot(TMMres, aes(x=Expected_LFC, y=logFC)) +
  annotate(geom = "text", y = 3, x = -3, label = eq_r2) +
  geom_point(aes(color = Spikes), shape = 1) +
  geom_smooth(method='lm', se=TRUE, color = 'blue', linetype = 'dashed') +
  theme_classic() +
  theme(aspect.ratio = 1) +
  xlab("Expected LFC") +
  ylab("Observed LFC") +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed')
ggsave("project_results/metasequin_concentration/MixA_vs_MixB_LFC_TMMnormalise_edgeR.pdf")

# Plot RLE with RUVg normalisation
spikes <- rownames(total_counts)[grepl("^ONT", rownames(total_counts))]
set1 <- RUVg(set, spikes, k=1)

pdf("project_results/metasequin_concentration/RLE_TMM_RUVg_normalised.pdf")
plotRLE(set1, outline=FALSE)
dev.off()

# Plot RLE with TMM and RUVg normalisation
y <- calcNormFactors(y, method="RLE" )
design <- model.matrix(~sampleType + pData(set1)$W_1)
y <- estimateDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
RUVgres <- data.frame(topTags(lrt, n = Inf))

# Add expected LFC information
idx <- match(rownames(RUVgres), sequin_conc$ID)
RUVgres$MIX_A <- sequin_conc$MIX_A [idx]
RUVgres$MIX_B <- sequin_conc$MIX_B [idx]
RUVgres$MIX_A[is.na(RUVgres$MIX_A)] <- 1
RUVgres$MIX_B[is.na(RUVgres$MIX_B)] <- 1
RUVgres$Expected_LFC <- log2(RUVgres$MIX_B/RUVgres$MIX_A)
RUVgres$Spikes <- grepl("^ONT", rownames(RUVgres))

# Linear model of observed vs expected
mod <- lm(logFC ~ Expected_LFC, RUVgres)
eq <- paste0("y = ", format(unname(coef(mod)[1]), digits = 2), " + ", format(unname(coef(mod)[2]), digits = 2), "x")
r2 <- paste0("r2 ", format(summary(mod)$r.squared, digits = 3))
eq_r2 <- paste0(eq, " ", r2)

# Plot observed vs expected
ggplot(RUVgres, aes(x=Expected_LFC, y=logFC)) +
  annotate(geom = "text", y = 3, x = -3, label = eq_r2) +
  geom_point(aes(color = Spikes), shape = 1) +
  geom_smooth(method='lm', se=TRUE, color = 'blue', linetype = 'dashed') +
  theme_classic() +
  theme(aspect.ratio = 1) +
  xlab("Expected LFC") +
  ylab("Observed LFC") +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed')
ggsave("project_results/metasequin_concentration/MixA_vs_MixB_LFC_RUVgnormalise_edgeR.pdf")

