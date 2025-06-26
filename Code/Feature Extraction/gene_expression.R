library(dplyr)
library(tidyverse)
library(readr)
library(data.table)
library(tidyr)

# Calculate measures of gene expression diversity (Tau scores and Shannon entropy) using BioQC
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BioQC")
library(BioQC)
??BioQC

# Working with high throughput gene expression data from FB_2025_02 release
# Has both FlyAtlas2 (FPKM) and modENCODE (RPKM) data for each gene 
# 2981421 obs. of 9 variables
expression_data <- fread("~/FYP/data/high-throughput_gene_expression_fb_2025_02.tsv", skip = 5)
EG_ids <- read.table('~/FYP/Flybase_IDs/essential_genes.txt')
NEG_ids <- read.table('~/FYP/Flybase_IDs/nonessential_genes.txt')

# Extract data for genes of interest and segregate FlyAtlas and modENCODE datasets
# Note, Knoblich Neural Cell RNA-seq, Lai_miRNA, and Casas-Vila datasets are not included in analysis
filtered_expression <- na.omit(expression_data[Gene_ID %in% NEG_ids$V1]) # 373793 x 9 
unique(expression_data$Dataset_Name) 

# Run first three only
flyatlas_data <- filtered_expression[Dataset_Name == "FlyAtlas2"] # 1247 x 9 
modencode_tissue <- filtered_expression[Dataset_Name == "modENCODE_mRNA-Seq_tissues"] # 841 x 9
modencode_development <- filtered_expression[Dataset_Name == "modENCODE_mRNA-Seq_development"] # 870 x9
modencode_cellb <- filtered_expression[Dataset_Name == "modENCODE_mRNA-Seq_cell.B"] # 696 x 9
modencode_treatments <- filtered_expression[Dataset_Name == "modENCODE_mRNA-Seq_treatments"] # 783 x 9

# Exploring FlyAtlas2 and ModEncode datasets
# Seems to be data for every condition for every gene (eg gene X has expression reported across 43 tissues)
unique(flyatlas_data$Sample_Name) # 43 different tissues (43*29 genes = 1247, which is the size of the df)
unique(modencode_tissue$Sample_Name) # 29 different tissues 
unique(modencode_development$Sample_Name) # 30 different developmental stages
unique(modencode_cellb$Sample_Name) # 24 different cell lines
unique(modencode_treatments$Sample_Name) # 27 different treatments 

# For my investigation, tissue and developmental-stage data seems most relevant 
# Cell line and treatment data could be useful for exploring conditional essentiality but 
# unnecessary for my FYP. 

# Organising FlyAtlas and modEncode
flyatlas_data <- flyatlas_data[, .(Sample_Name, Gene_ID, Gene_Symbol, Expression_Value)]
modencode_tissue <- modencode_tissue[, .(Sample_Name, Gene_ID, Gene_Symbol, Expression_Value)]
modencode_development <- modencode_development[, .(Sample_Name, Gene_ID, Gene_Symbol, Expression_Value)]

# Apply Log2 transform to FPKM and RPKM data
flyatlas_data$log_expression <- log2(flyatlas_data$Expression_Value + 1)
modencode_tissue$log_expression <- log2(modencode_tissue$Expression_Value + 1)
modencode_development$log_expression <- log2(modencode_development$Expression_Value + 1)

# Summary stats for flyatlas2 using log transformed data 
flyatlas_summary <- flyatlas_data[, .(
  fa_mean_expression = mean(log_expression, na.rm = TRUE),
  fa_median_expression = median(log_expression, na.rm = TRUE),
  fa_max_expression = max(log_expression, na.rm = TRUE),
  fa_tissue_with_max_expression = Sample_Name[which.max(log_expression)],
  fa_expression_variance = var(log_expression, na.rm = TRUE),
  fa_expression_sd = sd(log_expression, na.rm = TRUE), 
  fa_expression_breadth = sum(log_expression > 1, na.rm = TRUE)
), by = Gene_ID]

# Re-format data 
flyatlas_wide <- dcast(flyatlas_data, Gene_ID ~ Sample_Name, value.var = "log_expression", 
                       fun.aggregate = mean)
flyatlas_expr_matrix <- as.matrix(flyatlas_wide[, -1]) # Convert to matrix with rownames
rownames(flyatlas_expr_matrix) <- flyatlas_wide$Gene_ID

# Calculate entropy metrics from BioQC
flyatlas_entropy <- flyatlas_data %>%
  group_by(Gene_ID) %>%
  summarise(
    shannon_entropy = BioQC::entropy(log_expression)
  )

flyatlas_ES <- BioQC::entropySpecificity(flyatlas_expr_matrix, norm = FALSE) # [1:29]

# Define function to calculate tau score (metric of tissue-specific expression)
calculate_tau_score <- function(expr_matrix) {
  # Initialise an empty vector to store Tau scores
  tau_scores <- numeric(nrow(expr_matrix))
  
  # Loop over each gene (row in the matrix)
  for (i in 1:nrow(expr_matrix)) {
    gene_expr <- expr_matrix[i, ]
    
    # Skip genes that have no expression data (all zeros)
    if (all(gene_expr == 0)) {
      tau_scores[i] <- NA
    } else {
      # Calculate max expression, mean expression, and min non-zero expression
      max_expr <- max(gene_expr)
      mean_expr <- mean(gene_expr[gene_expr > 0])  # exclude zeros
      min_nonzero_expr <- min(gene_expr[gene_expr > 0])
      
      # Calculate Tau score using the formula
      tau_scores[i] <- (max_expr - mean_expr) / (max_expr - min_nonzero_expr)
    }
  }
  
  # Return a named vector of Tau scores
  names(tau_scores) <- rownames(expr_matrix)
  return(tau_scores)
}

# Apply function to calculate tau score
flyatlas_tau <- calculate_tau_score(flyatlas_expr_matrix)

# Organise into final FlyAtlas2 feature set (not including gene interactions)
flyatlas_features <- data.frame(
  Gene_ID = rownames(flyatlas_expr_matrix),
  fa_Tau_Score = flyatlas_tau,
  fa_Shannon_Entropy = flyatlas_entropy$shannon_entropy,
  fa_Shannon_Entropy_Specificity = flyatlas_ES
)

flyatlas_features <- right_join(flyatlas_summary, flyatlas_features, by = "Gene_ID")

# Repeat with modEncode tissue 
modencode_tissue_summary <- modencode_tissue[, .(
  met_mean_expression = mean(log_expression, na.rm = TRUE),
  met_median_expression = median(log_expression, na.rm = TRUE),
  met_max_expression = max(log_expression, na.rm = TRUE),
  met_tissue_with_max_expression = Sample_Name[which.max(log_expression)],
  met_expression_variance = var(log_expression, na.rm = TRUE),
  met_expression_sd = sd(log_expression, na.rm = TRUE), 
  met_expression_breadth = sum(log_expression > 1, na.rm = TRUE)
), by = Gene_ID]


# Organise expression matrix
met_wide <- dcast(modencode_tissue, Gene_ID ~ Sample_Name, value.var = "log_expression", 
                       fun.aggregate = mean)
met_expr_matrix <- as.matrix(met_wide[, -1]) # Convert to matrix with rownames
rownames(met_expr_matrix) <- met_wide$Gene_ID

# Calculate entropy metrics from BioQC
met_entropy <- modencode_tissue %>%
  group_by(Gene_ID) %>%
  summarise(
    shannon_entropy = BioQC::entropy(log_expression)
  )

# Calculate entropy specificity 
met_ES <- BioQC::entropySpecificity(met_expr_matrix, norm = FALSE) # [1:29]

# Calculate tau score 
met_tau <- calculate_tau_score(met_expr_matrix)

# Organise into modEncode tissue feature set 
modencode_tissue_features <- data.frame(
  Gene_ID = rownames(met_expr_matrix),
  met_Tau_Score = met_tau,
  met_Shannon_Entropy = met_entropy$shannon_entropy,
  met_Shannon_Entropy_Specificity = met_ES
)

modencode_tissue_features <- right_join(modencode_tissue_summary, 
                                        modencode_tissue_features, by = "Gene_ID")

# Repeat for modEncode Developmental Stage 
modencode_dev_summary <- modencode_development[, .(
  med_mean_expression = mean(log_expression, na.rm = TRUE),
  med_median_expression = median(log_expression, na.rm = TRUE),
  med_max_expression = max(log_expression, na.rm = TRUE),
  med_dev_stage_with_max_expression = Sample_Name[which.max(log_expression)],
  med_expression_variance = var(log_expression, na.rm = TRUE),
  med_expression_sd = sd(log_expression, na.rm = TRUE), 
  med_expression_breadth = sum(log_expression > 1, na.rm = TRUE)
), by = Gene_ID]

# Organise expression matrix 
med_wide <- dcast(modencode_development, Gene_ID ~ Sample_Name, value.var = "log_expression", 
                  fun.aggregate = mean)
med_expr_matrix <- as.matrix(med_wide[, -1]) # Convert to matrix with rownames
rownames(med_expr_matrix) <- med_wide$Gene_ID

# Calculate entropy metrics from BioQC
med_entropy <- modencode_development %>%
  group_by(Gene_ID) %>%
  summarise(
    shannon_entropy = BioQC::entropy(log_expression)
  )

# Calculate entropy specificity 
med_ES <- BioQC::entropySpecificity(med_expr_matrix, norm = FALSE) # [1:29]

# Calculate tau score 
med_tau <- calculate_tau_score(med_expr_matrix)

# Organise into modEncode tissue feature set 
modencode_dev_features <- data.frame(
  Gene_ID = rownames(med_expr_matrix),
  med_Tau_Score = med_tau,
  med_Shannon_Entropy = med_entropy$shannon_entropy,
  med_Shannon_Entropy_Specificity = med_ES
)

modencode_dev_features <- right_join(modencode_dev_summary, 
                                        modencode_dev_features, by = "Gene_ID")

# Save features 
write.csv(flyatlas_features, "~/FYP/feature_extraction/NEG_flyatlas_features.csv", row.names = FALSE)
write.csv(modencode_tissue_features, "~/FYP/feature_extraction/NEG_met_features.csv", row.names = FALSE)
write.csv(modencode_dev_features, "~/FYP/feature_extraction/NEG_med_features.csv", row.names = FALSE)
