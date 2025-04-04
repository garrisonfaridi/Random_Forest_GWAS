#!/bin/bash
#SBATCH --job-name=plink_pca      # Job name
#SBATCH --output=logs/plink_pca.out    # Output log file
#SBATCH --error=logs/plink_pca.err     # Error log file
#SBATCH --time=2:00:00            # Max runtime 
#SBATCH --nodes=1                 # Use a single node
#SBATCH --ntasks=1                # Single task
#SBATCH --cpus-per-task=4         # 4 CPU cores for processing
#SBATCH --mem=32G                 # Request 32GB RAM 

# Load PLINK module 
module load plink/1.90b6.21   # Adjust based on available PLINK version

# Define input and output filenames
INPUT_FILE="genotypes/genotypes_816k_bin/genotypes_816k"   # Input PLINK binary dataset (geno_binary.bed, .bim, .fam)
#FILTERED_FILE="genotypes/genotypes_maf0.05_miss0.1_bin/genotypes_maf0.03_miss0.1"  # Filtered genotypes file
PCA_OUTPUT="PCAs/pca_results_816k"   # Output filename for PCA results
FINAL_PCA_OUTPUT="${PCA_OUPUT}_10pcs.txt"

# Step 1: Filter SNPs and individuals with >10% missingness (filter if your binary genotypes file isn't already)
# plink --bfile $INPUT_FILE --maf 0.03 --geno 0.1 --make-bed --out $FILTERED_FILE

# Step 2: Perform PCA with 10 components on the cleaned dataset
plink --bfile $INPUT_FILE --pca 10 --out $PCA_OUTPUT 



# Extract the first 10 PCs for downstream analysis
echo -e "FID IID PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10" > $FINAL_PCA_OUTPUT
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' ${PCA_OUTPUT}.eigenvec >> ${FINAL_PCA_OUTPUT}

# Completion message
echo "PLINK missingness filtering and PCA completed successfully. Results saved as ${FINAL_PCA_OUTPUT}"