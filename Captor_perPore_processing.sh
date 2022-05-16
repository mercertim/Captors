# --------------------------------------------------------------------------
# Generate Captor error rate per pore
# --------------------------------------------------------------------------

# # Package requirements (command line)
# samtools
# bamtools
# 
# # Custom scripts
# script.py
# merge.py

# Paths
outDir = path/to/output
refDir = path/to/reference_sequences

# Make array of samples
sample_arr=(Sample_1 Sample_2 Sample_3)

# Loop across each sample
for sample in ${sample_arr[@]}; do

  # Split Captor aligned bam into adaptors
  bamtools split -in $outDir/$sample/$sample'_adaptor.bam' -reference
  
  # Move output to dedicated directory
  mkdir $outDir/$sample/split_bam
  mv $outDir/$sample/$sample'_adaptor.REF'* $outDir/$sample/split_bam/
  rm $outDir/$sample/$sample'_adaptor.REF_unmapped.bam'
  
  # Move the .bam file for each Captor to its own directory
  while read line; do
    mkdir $outDir/$sample/split_bam/$line;
    mv $outDir/$sample/split_bam/*$line.bam $outDir/$sample/split_bam/$line
  done < $refDir/Captor_names.txt
  
  # Sort and index .bam files for each Captor
  while read line; do
    samtools sort $outDir/$sample/split_bam/$line/*$line.bam > $outDir/$sample/split_bam/$line/$line'_sorted.bam'
    samtools index $outDir/$sample/split_bam/$line/$line'_sorted.bam'
  done < $refDir/Captor_names.txt
  
  # Make single line bed files for each Captor to subset .bams
  while read line; do
    grep $line $refDir/Captors.bed > $outDir/$sample/split_bam/$line/$line'_one_line_bed.bed'
  done < $refDir/Captor_names.txt
  
  # Generate per pore stats for each Captor
  while read line; do
    mkdir $outDir/$sample/split_bam/$line/$line'_perRead'
    
    samtools view -b -L $outDir/$sample/split_bam/$line/$line'_one_line_bed'.bed $outDir/$sample/split_bam/$line/$line'_sorted.bam' > $outDir/$sample/split_bam/$line/$line'_perRead'/$line'_sorted_one_line_intersect.bam'
    samtools index $outDir/$sample/split_bam/$line/$line'_perRead'/$line'_sorted_one_line_intersect.bam'
    
    # Generate per base error stats for each read
    python3 script.py $refDir/Captors.fa $outDir/$sample/$sample'_combined_first_500.fastq' $outDir/$sample/split_bam/$line/$line'_one_line_bed'.bed $outDir/$sample/split_bam/$line/$line'_perRead'/$line'_sorted_one_line_intersect.bam'  > $outDir/$sample/split_bam/$line/$line'_perRead'/$line'_sorted_one_line_intersect'.perBase.tsv
    
    # Generate error stats per pore
    python3 merge.py $outDir/$sample/split_bam/$line/$line'_perRead'/$line'_sorted_one_line_intersect'.perBase.tsv > $outDir/$sample/split_bam/$line/$line'_perRead'/$line'_sorted_one_line_intersect'.perRead.tsv
  done < $refDir/Captor_names.txt
  
done
