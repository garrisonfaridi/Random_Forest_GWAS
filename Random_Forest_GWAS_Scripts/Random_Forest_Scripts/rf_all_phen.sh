#!/bin/bash
#SBATCH --job-name=rf_all_phen   # Job name
#SBATCH --output=logs/rf_all_phen.out  # Output log file
#SBATCH --error=logs/rf_all_phen.err   # Error log file
#SBATCH --time=12:00:00   # Max runtime (adjust as needed)
#SBATCH --nodes=1         # Use a single node
#SBATCH --ntasks=1        # Single task
#SBATCH --cpus-per-task=12  # 4 CPU cores 
#SBATCH --mem=24G         # 24GB RAM 
#SBATCH --mail-user=gpf2024@nyu.edu  # Replace with your email
#SBATCH --mail-type=BEGIN,END,FAIL

# Load R module if using an HPC system
module load r/gcc/4.4.0 

# Run the R script
Rscript R_scripts/rf_gwas_allPhen.R
