#!/bin/bash

#SBATCH --job-name=alignment # you can give your job a name
#SBATCH --ntasks=24 # the number of processors or tasks
#SBATCH --account=itcga # our account
#SBATCH --reservation=ITCGA2025 # this gives us special access during the workshop
#SBATCH --time=10:00:00 # the maximum time for the job
#SBATCH --mem=32gb # the amount of RAM
#SBATCH --partition=itcga # the specific server in chimera we are using
#SBATCH --error=log/%x-%A.err   # a filename to save error messages into
#SBATCH --output=log/%x-%A.out  # a filename to save any printed output into

module load gcc-10.2.0-gcc-9.3.0-f3oaqv7
module load python-3.8.12-gcc-10.2.0-oe4tgov
module load hisat2-2.1.0-gcc-9.3.0-u7zbyow

# Define variables
index_dir=$1
input_dir=$2 # takes this from the command line, first item after the script
output_dir=$3 # takes this from the command line, second item

# Loop through R1 files in the input directory
for file in "$input_dir"/*_1_trim.fastq.gz; 
do
  # Extract the base name without the suffix (e.g., "C1_S4_L001")
  base=$(basename "$file" _1_trim.fastq.gz)

  hisat2 -p 24 \
  -x "$index_dir" \
  -1 "$input_dir/${base}_1_trim.fastq.gz" \
  -2 "$input_dir/${base}_2_trim.fastq.gz" \
  -S "$output_dir/${base}.sam"

done

echo "Finish Run"
echo "End time: $(date)"
