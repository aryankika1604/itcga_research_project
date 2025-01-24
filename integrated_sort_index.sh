#!/bin/bash
#SBATCH --job-name=sort_and_index_bam
#SBATCH --ntasks=24
#SBATCH --account=itcga
#SBATCH --reservation=ITCGA2025
#SBATCH --time=10:00:00
#SBATCH --mem=32gb
#SBATCH --partition=itcga
#SBATCH --error=log/%x-%A.err
#SBATCH --output=log/%x-%A.out

module load samtools-1.10-gcc-9.3.0-flukja5

# Define the input directory
input_dir=$1

# Loop through each BAM file in the input directory
for file in "$input_dir"/*.bam; do
    # Extract the base name without the suffix
    base=$(basename "$file" .bam)

    # Define temporary output for sorted BAM
    sorted_bam="$input_dir/${base}_sorted.bam"

    # Sort the BAM file and write the sorted output to a pipe
    # Simultaneously index the BAM file
    echo "Processing file: $file"
    samtools sort "$file" -o "$sorted_bam" && \
    samtools index "$sorted_bam"

    if [[ $? -eq 0 ]]; then
        echo "Successfully processed $file"
    else
        echo "Error processing $file"
        exit 1
    fi
done

echo "Finish Run"
echo "End time: $(date)"

