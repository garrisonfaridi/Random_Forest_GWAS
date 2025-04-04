# RF GWAS 400 trees
# Load necessary libraries
library(randomForest)
library(data.table)

# Load genotype data
geno <- fread("genotypes/genotypes_816k_bin/genotypes_816k.raw.raw")


# Write Mode Imputation function
replace_na_with_mode <- function(dt) {
  for (j in seq_along(dt)) {
    col <- dt[[j]]
    if (is.numeric(col)) {
      mode_value <- as.numeric(names(sort(table(col), decreasing = TRUE)[1]))  # Fast mode calculation
      set(dt, which(is.na(col)), j, mode_value)  # Fast in-place NA replacement
    }
  }
}

# Apply mode imputation directly to data.table
replace_na_with_mode(geno)

pca <- fread("PCAs/816k_10pcs.txt")  # PCA components
pheno <- fread("Phenotypes/dat_tot_ForGemmaGF_300.csv")  # Phenotype data

# Merge All Data
# Add FID column to phenotype data to match genotypes
pheno[, FID := geno$FID]
setcolorder(pheno, c("FID", setdiff(names(pheno), "FID")))

# Remove rows that are completely NA
pheno <- pheno[rowSums(is.na(pheno[, 3:ncol(pheno), with = FALSE])) < (ncol(pheno) - 2)]

# Merge all at once
merged_data <- merge(merge(geno, pca, by = c("FID", "IID")), pheno, by = "FID")

# Remove metadata columns and phenotypes while merging (faster than subsetting later)
phenotype_cols <- colnames(pheno)[-1]  # Exclude 'FID'
X_816k <- merged_data[, !c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE", phenotype_cols), with = FALSE]


# Create empty list to store results
all_importances <- vector("list", length(colnames(pheno)[-1]))

# Precompute valid rows for each phenotype
valid_rows_list <- lapply(colnames(pheno)[-1], function(phenotype_col) {
  complete.cases(merged_data[[phenotype_col]])
})

# Loop over all 17 phenotype columns
for (i in seq_along(colnames(pheno)[-1])) {
  
  phenotype_col <- colnames(pheno)[-1][i]  # Exclude 'FID'
  
  # Pre-filter valid data (computed earlier)
  valid_rows <- valid_rows_list[[i]]
  X <- X_816k[valid_rows, ]
  y <- merged_data[[phenotype_col]][valid_rows]
  
  # Train Random Forest Model
  rf_model <- randomForest(x = X, y = y, ntree = 400, mtry = floor(sqrt(ncol(X))), importance = TRUE)
  
  # Extract feature importance
  importance_df <- as.data.frame(importance(rf_model))
  importance_df$Feature <- rownames(importance_df)
  importance_df$Phenotype <- phenotype_col  # Add phenotype name
  
  # Save importance data
  all_importances[[i]] <- importance_df
  
  print(paste("RF for", phenotype_col, "completed"))
}

# Combine all results into a single dataframe
final_importance_df <- rbindlist(all_importances)

# Save feature importance to a CSV file
fwrite(final_importance_df, "feature_importance_816k_all_phenotypes_SLURM.csv")

# Print completion message
cat("Random Forest models completed for all 17 phenotypes.\nFeature importance saved to 'feature_importance_816k_all_phenotypes.csv'.\n")
