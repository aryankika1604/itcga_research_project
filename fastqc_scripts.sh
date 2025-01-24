#!/bin/bash

#SBATCH --job-name=fastqc_Aryan # you can give your job a name
#SBATCH --nodes 1 # the number of processors or tasks
#SBATCH --cpus-per-task=2
#SBATCH --account=itcga # our account
#SBATCH --time=1:00:00 # the maximum time for the job
#SBATCH --mem=4gb # the amount of RAM
#SBATCH --partition=itcga # the specific server in chimera we are using
#SBATCH --reservation=ITCGA2025
#SBATCH --error=logs/%x-%A_%a.err   # error messages file with date and time
#SBATCH --output=logs/%x-%A_%a.out  # output file with date and time

# Loading the fastqc script from 'module'
module load fastqc-0.11.9-gcc-10.2.0-osi6pqc 

# Running the fastqc script on the files with .fastq extention in 1_project/data
fastqc -o /itcgastorage/data01/itcga_workshops/aug2024_genomics/team_cherries/results /itcgastorage/data01/itcga_workshops/aug2024_genomics/team_cherries/raw_data/*.fastq.gz

#mv /itcgastorage/share_home/a.kikaganeshwala001/1_project/data/*fastqc* /itcgastorage/share_home/a.kikaganeshwala001/1_project/results

echo "fastqc implemented successfully"

