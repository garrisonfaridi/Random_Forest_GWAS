# -------------------------------------------
# EXTRACT TOP 40 IMPORTANCES FOR EACH PHENOTYPE
# -------------------------------------------

library(data.table)
library(dplyr)


# Read feature importance results
importances_allPhen <- fread("feature_importance_816k_all_phenotypes_SLURM.csv")

# Separate into 17 tables by phenotype
# Define ouput directory
output_dir <- "Importance_Datasets"

# Keep only the the column headers of the phenotye key as a list
phenotype_list <- c(
  "7mSh",
  "8MTO",
  "3mSOp",
  "4mSOb",
  "5mSOp",
  "6mSOh",
  "7mSOh",
  "8mSOo",
  "Pren",
  "Buen",
  "Peen",
  "S2hBuen",
  "2hPeen",
  "IM",
  "1moIM",
  "1hIM",
  "4moIM"
)

# Split datset by Phenotype and save each sorted table
unique_phenotypes <- unique(importances_allPhen$Phenotype)

# Replace V1, V2, ... with phenotype_list values
if (length(unique_phenotypes) == length(phenotype_list)) {
  importances_allPhen$Phenotype <- phenotype_list[match(importances_allPhen$Phenotype, unique_phenotypes)]
} else {
  stop("Error: The number of unique phenotypes in the dataset does not match the length of phenotype_list.")
}

# Split the dataset by phenotype and save each osrted table

for (phenotype in unique(phenotype_list)) {
  phenotype_data <- importances_allPhen %>% 
    filter(Phenotype == phenotype) %>%
    arrange(desc(IncNodePurity))
  
  # Write csv
  write_csv(phenotype_data, file.path(output_dir, paste0("importance_", phenotype, ".csv"))
  )
  
  # Print message
  print(paste("Saved importance data for", phenotype))
}

# Save Top 40 Importances for Each Phenotype
input_dir <- "Importance_Datasets"
output_dir <- "Importance_Top40"

# Create output dir
dir.create(output_dir, showWarnings = FALSE)

# Regular expression function to extract `ps`
extract_ps <- function(feature) {
  if (grepl("^PC[0-9]+$", feature)) {
    return(NA)  # If the feature is PC1, PC2, etc., return NA
  } else {
    feature <- gsub("^[0-9]+:", "", feature)  # Remove chromosome prefix (e.g., "1:", "2:")
    feature <- gsub("_[ATCG]$", "", feature)  # Remove allele suffix (e.g., "_A", "_T", "_G", "_C")
    return(feature)
  }
}
# Loop through each phenotype, read the dataset, subset top 40%, and save the file
for (phenotype in phenotype_list) {
  file_path <- file.path(input_dir, paste0("importance_",phenotype,".csv"))
  
  # Read the dataset
  if (file.exists(file_path)) {
    phenotype_data <- read_csv(file_path)
    
    # Create new "ps" column extracting position value
    phenotype_data <- phenotype_data %>%
      mutate(ps = sapply(Feature, extract_ps))
    
    # Determine the number ofrows corresponding to the top 40%
    top_n <- ceiling(0.40 * nrow(phenotype_data))
    
    # Supbet the top 40% based on IncNodePurity
    top_40_data <- phenotype_data %>%
      arrange(desc(IncNodePurity)) %>%
      head(top_n)
    
    # Save the subset dataset to output dir
    write.csv(top_40_data, file.path(output_dir, paste0("importance_", phenotype, "_top40.csv")))
  }
}