#!/bin/bash
#SBATCH --ntasks=2
#
#SBATCH --job-name=trim # you can give your job a name
#SBATCH --nodes 1 # the number of processors or tasks
#SBATCH --cpus-per-task=4
#SBATCH --account=itcga # our account
#SBATCH --time=4:00:00 # the maximum time for the job
#SBATCH --mem=4gb # the amount of RAM
#SBATCH --partition=itcga # the specific server in chimera we are using
#SBATCH --reservation=ITCGA2025
#SBATCH --error=logs/%x-%A_%a.err   # error messages file with date and time
#SBATCH --output=logs/%x-%A_%a.out  # output file with date and time

module load py-cutadapt-2.10-gcc-10.2.0-2x2ytr5
module load py-dnaio-0.4.2-gcc-10.2.0-gaqzhv4
module load py-xopen-1.1.0-gcc-10.2.0-5kpnvqq

#cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --trim-n -m 25
#-o /itcgastorage/share_home/a.kikaganeshwala001/1_project/data/C1_S4_L001_R1_001_downsampled_trimmed.fastq
#-p /itcgastorage/share_home/a.kikaganeshwala001/1_project/data/C1_S4_L001_R2_001_downsampled_trimmed.fastq
#/itcgastorage/share_home/a.kikaganeshwala001/1_project/data/C1_S4_L001_R1_001_downsampled.fastq
#/itcgastorage/share_home/a.kikaganeshwala001/1_project/data/C1_S4_L001_R2_001_downsampled.fastq

# Define variables.
input_dir=$1 # takes this from the command line, first item after the script
output_dir=$2 # takes this from the command line, second item

for file in ${input_dir}/*.fastq.gz
 do
 # Pull basename
 name=$(basename ${file} _1.fastq.gz)

 # Run cutadapt
 cutadapt --cores=2  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
 -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --trim-n -m 25 \
 -o ${output_dir}/${name}_1_trim.fastq.gz \
 -p ${output_dir}/${name}_2_trim.fastq.gz \
 ${input_dir}/${name}_1.fastq.gz \
 ${input_dir}/${name}_2.fastq.gz

 echo cutadapt is finished with $name

done
