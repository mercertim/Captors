# --------------------------------------------------------------------------
# Generate per base error statistics for Captors and Metasequins
# --------------------------------------------------------------------------

# # Package requirements (command line)
# fastp
# minimap2
# samtools
# 
# # Package requirements (python 3)
# pysamstats
# 
# # Custom scripts
# analyzePile.py

# Paths
inDir=. #path/to/sample
outDir=. #path/to/output
refDir=. #path/to/reference_sequences

# Get per base error statistics for Captors
# --------------------------------------------------------------------------

# Create sample array
sample_arr=(demo_data)

# Loop for all samples in array
for sample in ${sample_arr[@]}; do

  # Make output directory
  mkdir $outDir/$sample
  
  # Trim reads to first 500bp
  fastp --disable_adapter_trimming --disable_trim_poly_g --disable_quality_filtering --max_len1 500 --thread 16 -i $inDir/$sample.fastq -o $outDir/$sample/$sample'_combined_first_500.fastq'
  
  # Align to Captor sequence with minimap2 and sort bam
  minimap2 -ax map-ont -t 8 $refDir/Captors.fa $outDir/$sample/*'_combined_first_500.fastq' | samtools sort - > $outDir/$sample/$sample'_captor.bam'
  
  # Index bam file
  samtools index $outDir/$sample/$sample'_captor.bam'
  
  # Get pileup stats per base for each adaptor
  pysamstats --fasta $refDir/Captors.fa --type variation $outDir/$sample/$sample'_captor.bam' > $outDir/$sample/$sample'_captor.bam.bed'
  
  # Collate pileup stats in dataframe
  python3 analyzePile.py $outDir/$sample/$sample'_captor.bam.bed' > $outDir/$sample/$sample'_captor.bam.bed.tsv'
  
done

# Get per base error statistics for Metasequins (replace with sample genome)
# --------------------------------------------------------------------------

# Loop across the same samples from above array (same as Captors)
for sample in ${sample_arr[@]}; do

  # Align to Metasequins sequence with minimap2 and sort bam
  minimap2 -ax map-ont -t 8 $refDir/MetaSequin.fa $inDir/$sample.fastq | samtools sort - > $outDir/$sample/$sample'_metasequin.bam'
  
  # Index bam file
  samtools index $outDir/$sample/$sample'_metasequin.bam'
  
  # Get pileup stats per base for each Metasequins
  pysamstats --fasta $refDir/MetaSequin.fa --type variation $outDir/$sample/$sample'_metasequin.bam' > $outDir/$sample/$sample'_metasequin.bam.bed'
  
  # Collate pileup stats in dataframe
  python3 analyzePile.py $outDir/$sample/$sample'_metasequin.bam.bed' > $outDir/$sample/$sample'_metasequin.bam.bed.tsv'
  
  # Finish loop
done
