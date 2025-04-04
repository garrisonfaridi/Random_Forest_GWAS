#!/bin/bash
#SBATCH --job-name=plink_bin2raw
#SBATCH --output=logs/plink_bin2raw_output.log
#SBATCH --error=logs/plink_bin2raw_error.log
#SBATCH --time=72:00:00   # Set max job time
#SBATCH --mem=12G  # Adjust memory as needed
#SBATCH --cpus-per-task=4 #CPU cores
#SBATCH --ntasks=1  # Number of cores/tasks

# Load Plink 
module load plink/1.90b6.21 

# Convert binary file to raw
plink --bfile genotypes/genotypes_816k_bin/genotypes_816k --recodeA --out genotypes/genotypes_816k_bin/genotypes_816k.raw