#!/bin/bash
#SBATCH --job-name=plink_extract
#SBATCH --output=logs/plink_extract_output.log
#SBATCH --error=logs/plink_extract_error.log
#SBATCH --time=72:00:00   # Set max job time
#SBATCH --mem=10G  # Adjust memory as needed
#SBATCH --cpus-per-task=1 #CPU cores
#SBATCH --ntasks=1  # Number of cores/tasks

module load plink/1.90b6.21 

plink --bfile genotypes/genotypes_305_bin/genotypes_305 --extract genotypes/snp_list_816k.txt --make-bed --out genotypes_816k

