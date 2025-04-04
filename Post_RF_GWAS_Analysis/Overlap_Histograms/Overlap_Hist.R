# -------------------------------------------------------------------
# Make histograms of overlapping SNP counts over genes
# -------------------------------------------------------------------

# Load necessary library
library(ggplot2)

# Read the dataset
ali <- read.csv("aliphatic_snp_overlap_avg.csv", stringsAsFactors = FALSE)  # Replace with your actual file
ind <- read.csv("indolic_snp_overlap_avg.csv", stringsAsFactors = FALSE)  # Replace with your actual file

# Create a bar plot
ggplot(ali, aes(x = gene_abbreviation, y = num_snps_avg)) +
  geom_bar(stat = "identity", fill = "steelblue") +  # Use bar plot with actual values
  theme_minimal() +
  labs(title = "Average SNP Overlap per Gene Aliphatic",
       x = "Gene Abbreviations",
       y = "Average Number of SNPs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability

# Create a bar plot
ggplot(ind, aes(x = gene_abbreviation, y = num_snps_avg)) +
  geom_bar(stat = "identity", fill = "steelblue") +  # Use bar plot with actual values
  theme_minimal() +
  labs(title = "Average SNP Overlap per Gene Indolic",
       x = "Gene Abbreviations",
       y = "Average Number of SNPs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability