#!/bin/bash

#SBATCH --job-name=gene_count # you can give your job a name
#SBATCH --ntasks=24 # the number of processors or tasks
#SBATCH --account=itcga # our account
#SBATCH --reservation=ITCGA2025 # this gives us special access during the workshop
#SBATCH --time=10:00:00 # the maximum time for the job
#SBATCH --mem=32gb # the amount of RAM
#SBATCH --partition=itcga # the specific server in chimera we are using
#SBATCH --reservation=ITCGA2025
#SBATCH --error=log/%x-%A.err   # a filename to save error messages into
#SBATCH --output=log/%x-%A.out  # a filename to save any printed output into

module load subread-2.0.2-gcc-10.2.0

# Define variables
annotation_dir=$1
output_dir=$2 # takes this from the command line, first item after the script
input_dir=$3 # takes this from the command line, second item

# Loop through R1 files in the input directory
for file in "$input_dir"/*_sorted.bam; do
  # Extract the base name without the suffix (e.g., "C1_S4_L001")
  base=$(basename "$file" _sorted.bam)

featureCounts -a "$annotation_dir/Homo_sapiens.GRCh38.111.gtf" \
 -o "$output_dir/${base}_counts.txt" \
 -T 24 \
 -p -B -C "$file"

done
