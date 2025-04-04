# ------------------------------------------------
# Calculate Candidate Gene SNP Overlap with Top 40 Percent Importance Scores for all phenotypes
# ------------------------------------------------

library(data.table)

##### DEFINE SNP LIST, CANDIDATE GENES #####

# load SNPs, maf > 0.05, TOU-A
snps.toua        = subset(read.delim("300_gluc1_maf0.03_miss0.1GF.assoc.txt", sep="\t")[,c("chr","ps","af")], af > 0.05)

# get full list of SNPs
snps.all         = snps.toua[,c("chr","ps")]

# rs from vcf sometimes has chr first (tou), sometimes ps first (regmap/1001 imputed), so standardize
snps.all$rs      = paste(snps.all$chr, snps.all$ps, sep = "_")
snps.all         = unique(snps.all)

# # save SNPs list so it can be loaded directly later
# write.csv(snps.all, "full_SNP_list_maf05_wTouFrachonMiss10.csv")
# # to load:
# snps.all = read.csv("full_SNP_list_maf05_wTouFrachonMiss10.csv")

# load candidate gene list
candidate.genes  = read.delim("pathway_BiosyntheticGenes_HarunTable2.txt", sep = "\t")
candidate.genes  = subset(candidate.genes, aliphatic == "yes" | indolic == "yes")
nrow(candidate.genes) # 45 candidate genes

# load table of gene models (left and right boundaries)
gene.models      = read.delim("TAIR10_GFF3_geneBoundariesTableADG_extendedWindows.txt", sep = "\t")

# merge gene model info with candidate gene list
candidate.genes  = merge(candidate.genes, gene.models, by = "gene_id")
nrow(candidate.genes) # 45 genes still retained with gene model coordinates now included
# head(candidate.genes, 10) # check gene model coordinates (start and end positions in bp on chromosome) for candidate genes

##### ALIPHATICS #####

candidate.genes.ali = subset(candidate.genes, aliphatic == "yes")

candidate.genes.ind = subset(candidate.genes, indolic == "yes")

candidate.snps.per.gene = list()

# get list of candidate SNPs per gene
for (i in 1:nrow(candidate.genes.ali)){
  
  snps = subset(snps.all, chr == candidate.genes.ali[i,"chr"] &
                  ps > candidate.genes.ali[i,"ps_left"] &
                  ps < candidate.genes.ali[i,"ps_right"],
                select = rs
  )
  
  candidate.snps.per.gene[[ as.character(candidate.genes.ali[i,"gene_id"]) ]] = snps
  
}
# Aliphatic Gluc Implortances file paths List
ali_importances_paths <- list("Importance_Top40/importance_7mSh_top40.csv",
                              "Importance_Top40/importance_8MTO_top40.csv",
                              "Importance_Top40/importance_3mSOp_top40.csv",
                              "Importance_Top40/importance_4mSOb_top40.csv",
                              "Importance_Top40/importance_5mSOp_top40.csv",
                              "Importance_Top40/importance_6mSOh_top40.csv",
                              "Importance_Top40/importance_7mSOh_top40.csv",
                              "Importance_Top40/importance_8mSOo_top40.csv",
                              "Importance_Top40/importance_Pren_top40.csv",
                              "Importance_Top40/importance_1hIM_top40.csv",
                              "Importance_Top40/importance_Buen_top40.csv",
                              "Importance_Top40/importance_Peen_top40.csv",
                              "Importance_Top40/importance_IM_top40.csv")

# Indolic Gluc Importances file paths list
ind_importances_paths <- list("Importance_Top40/importance_IM_top40.csv",
                              "Importance_Top40/importance_1moIM_top40.csv",
                              "Importance_Top40/importance_1hIM_top40.csv",
                              "Importance_Top40/importance_4moIM_top40.csv")

# Function to count overlapping SNPs per gene
count_snps_per_gene <- function(importances_file, candidate_genes_subset) {
  importances <- fread(importances_file)
  importances$ps <- as.numeric(importances$ps)  # Ensure `ps` column is numeric
  
  # Remove NAs from ps before filtering
  importances <- importances[!is.na(importances$ps), ]
  
  # Initialize empty dataframe for results
  snp_counts <- data.frame(gene_id = character(),
                           gene_abbreviation = character(),
                           num_snps = integer(),
                           stringsAsFactors = FALSE)
  
  # Loop through each candidate gene
  for (i in 1:nrow(candidate_genes_subset)) {
    gene <- candidate_genes_subset[i, ]
    
    # Skip iteration if gene boundaries are missing
    if (is.na(gene$ps_left) | is.na(gene$ps_right)) {
      next
    }
    
    # Find SNPs within gene boundaries (Fixed incorrect pipe usage)
    snps_in_gene <- importances %>%
      filter(ps >= gene$ps_left & ps <= gene$ps_right)
    
    # Ensure `snps_in_gene` exists, even if empty
    num_snps <- ifelse(nrow(snps_in_gene) > 0, nrow(snps_in_gene), 0)
    
    # Store the count
    snp_counts <- rbind(snp_counts, data.frame(
      gene_id = gene$gene_id,
      gene_abbreviation = gene$gene_abbreviation,
      num_snps = num_snps  # Use computed num_snps
    ))
  }
  
  return(snp_counts)
}

# Define lists for output storage
ali_results <- list()
ind_results <- list()

# Process aliphatic datasets
for (file_path in ali_importances_paths) {
  ali_results[[file_path]] <- count_snps_per_gene(file_path, candidate.genes.ali)
}

# Process indolic datasets
for (file_path in ind_importances_paths) {
  ind_results[[file_path]] <- count_snps_per_gene(file_path, candidate.genes.ind)
}

# Save results to CSV files
for (file_path in names(ali_results)) {
  fwrite(ali_results[[file_path]], paste0("ali_snps_overlap", basename(file_path)))
}
for (file_path in names(ind_results)) {
  fwrite(ind_results[[file_path]], paste0("ind_snps_overlap", basename(file_path)))
}