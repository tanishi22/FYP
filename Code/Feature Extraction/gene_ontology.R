# Extracting gene-level features with biomaRt 

install.packages("biomaRt")
library(biomaRt)
library(dplyr)
library(tidyr)

# Obtain Ensembl Genomes for D. mel, linked to FlyBase DB 
ensembl <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")

#If ensembl server is down, try alternative:
ensembl <- useEnsembl(
  biomart = "ensembl",
  dataset = "dmelanogaster_gene_ensembl"
)

# Load sterile_dom FlyBase IDs of protein-coding genes (n = 29)
sterile_dom_ids <- readLines("~/FYP/FlyBase_IDs/sterile_dom_pep_IDs.txt") 
sterile_rec <- read.csv("~/FYP/feature_extraction/SR/sterile_rec_isoforms.csv")
lethal_rec <- read.csv("~/FYP/feature_extraction/LR/lethal_rec_isoforms.csv")
neg <- read.csv("~/FYP/feature_extraction/NEG/neg_isoforms.csv")

sterile_rec_ids <- sterile_rec$FBgn_ID
lethal_rec_ids <- lethal_rec$FBgn_ID
neg_ids <- neg$FBgn_ID
gene_isoform_map


# Target gene and sequence ontology terms, prioritised from total list
ontology_attributes <- c("ensembl_gene_id", "ensembl_peptide_id" , 
                         "namespace_1003","goslim_goa_description")

EG_ids <- read.table("~/FYP/FlyBase_IDs/essential_genes.txt")
NEG_ids <- read.table("~/FYP/FlyBase_IDs/nonessential_genes.txt")

# Extract ontology terms # 8094 x 7 in raw table 
ontology_result <- getBM(
  attributes = ontology_attributes,
  filters = "flybase_gene_id",
  values = NEG_4_gene_isoform$FBgn_ID,
  mart = ensembl
) 
head(ontology_result)

# Split into batches of 500 for larger gene sets! 
batch_size <- 500
batches <- split(neg_ids, ceiling(seq_along(neg_ids) / batch_size))

ontology_result <- do.call(rbind, lapply(batches, function(batch) {
  getBM(
    attributes = ontology_attributes,
    filters = "flybase_gene_id",
    values = batch,
    mart = ensembl
  )
}))

length(unique(ontology_result$ensembl_gene_id))


# only extract ontology features and results for representative isoform 
# 5316 x 7 variables 
filtered_ontology_NEG4 <- ontology_result %>%
  filter(ensembl_peptide_id %in% NEG_4_gene_isoform$ensembl_peptide_id)

length(unique(filtered_ontology_NEG4$ensembl_gene_id))

EG_1_ontology <- filtered_ontology
EG_2_ontology <- filtered_ontology_EG2
NEG_1_ontology <- filtered_ontology_NEG1
NEG_2_ontology <- filtered_ontology_NEG2
NEG_3_ontology <- filtered_ontology_NEG3
NEG_4_ontology <- filtered_ontology_NEG4


## Filter to drop peptide ids 
EG_1_ontology$ensembl_peptide_id <- NULL
EG_2_ontology$ensembl_peptide_id <- NULL
NEG_1_ontology$ensembl_peptide_id <- NULL
NEG_2_ontology$ensembl_peptide_id <- NULL
NEG_3_ontology$ensembl_peptide_id <- NULL
NEG_4_ontology$ensembl_peptide_id <- NULL

# Feature 1: Count occurrences of each namespace per gene
counts_per_gene <- NEG_4_ontology %>%
  filter(namespace_1003 != "" & !is.na(namespace_1003)) %>%  # remove empty and NA
  group_by(ensembl_gene_id, namespace_1003) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = namespace_1003, values_from = count, values_fill = 0)

print(counts_per_gene)

# Feature 2: GOSlim terms as categorical featurees 
gene_terms <- NEG_4_ontology %>%
  filter(goslim_goa_description != "" & !is.na(goslim_goa_description)) %>%
  group_by(ensembl_gene_id) %>%
  summarise(goslim_terms = list(unique(goslim_goa_description)), .groups = 'drop')

list(unique(gene_terms))

# Concatenate into final ontology result for the gene set 
NEG_4_onto <- merge(counts_per_gene, gene_terms, by = 'ensembl_gene_id')
NEG_4_onto <- NEG_4_onto %>%
  mutate(goslim_terms = sapply(goslim_terms, function(x) paste(x, collapse = "; ")))

# Rename first column to gene ID for consistency
colnames(NEG_4_onto)[colnames(NEG_4_onto) == "ensembl_gene_id"] <- "FBgn_ID"

# save output
write.csv(NEG_4_onto, "~/FYP/feature_extraction/ontology/NEG_4_ontology.csv", row.names = FALSE)


# --- Process ontology features based on GO Slim. Trying to create a gene-by-GO Slim binary matrix
# Extract unique gene IDs and GO Slim terms
genes <- unique(filtered_ontology$ensembl_gene_id)
goslims <- unique(NEG_4_ontology$goslim_goa_description) # n = 59
length(goslims)

# make matrix with goslims
ontology_mat <- matrix(0, nrow = length(genes), ncol = length(goslims),
                       dimnames = list(genes, goslims))

# Fill in 1 where gene has GO Slim term
for(i in seq_len(nrow(filtered_ontology))) {
  gene <- filtered_ontology$ensembl_gene_id[i]
  goslim <- filtered_ontology$goslim_goa_description[i]
  ontology_mat[gene, goslim] <- 1
}

# Convert matrix to data frame and organise
goslim_df <- as.data.frame(ontology_mat)
goslim_df$ensembl_gene_id <- rownames(goslim_df) # Add gene IDs as a column
rownames(goslim_df) <- NULL
goslim_df <- goslim_df[, c(ncol(goslim_df), 1:(ncol(goslim_df) - 1))] # Reorder to put gene ID as first column
colnames(goslim_df) <- gsub(" ", "_", colnames(goslim_df)) #Rename columns
colnames(goslim_df)[-1] <- paste0("GOSlim_", colnames(goslim_df)[-1]) #add 'GOSlim_' prefix to all columns

# repeat for go_names 
# Get unique genes and GO names for matrix dimensions
go_names <- unique(filtered_ontology$name_1006) # n = 422 which isn't too bad tbh
length(go_names)

# Initialize zero matrix
go_ontology_mat <- matrix(0, nrow = length(genes), ncol = length(go_names),
                          dimnames = list(genes, go_names))

# Fill matrix, skipping NA safely
for (i in seq_len(nrow(filtered_ontology))) {
  gene <- filtered_ontology$ensembl_gene_id[i]
  go_name <- filtered_ontology$name_1006[i]
  
  # Skip if either gene or go_name is NA or empty
  if (is.na(gene) || is.na(go_name) || gene == "" || go_name == "") next
  
  # Trim whitespace so col names don't mess up matrix
  gene <- trimws(gene)
  go_name <- trimws(go_name)
  
  # Assign 1 if gene and GO name exist in matrix dims
  if (gene %in% rownames(go_ontology_mat) && go_name %in% colnames(go_ontology_mat)) {
    go_ontology_mat[gene, go_name] <- 1
  }
}

# Convert matrix to data frame and organise
go_df <- as.data.frame(go_ontology_mat)
go_df$ensembl_gene_id <- rownames(go_df) # Add gene IDs as a column
rownames(go_df) <- NULL
go_df <- go_df[, c(ncol(go_df), 1:(ncol(go_df) - 1))] # Reorder to put gene ID as first column
colnames(go_df) <- gsub(" ", "_", colnames(go_df)) #Rename columns
colnames(go_df)[-1] <- paste0("GO_", colnames(go_df)[-1]) #add 'GO_' prefix to all columns

binary_go <- merge(go_df, goslim_df, by = 'ensembl_gene_id') # 482 attributes x 29 genes
summary(binary_go)
write.csv(binary_go, "~/FYP/feature_extraction/NEG/NEG_binary_GO.csv", 
          row.names = TRUE)

# comparison between phenotypes
binary_go_LD <- read.csv("/Users/tanishi/FYP/feature_extraction/LD/LD_binary_GO.csv")
binary_go_SD <- read.csv("/Users/tanishi/FYP/feature_extraction/SD/SD_binary_GO.csv")
binary_go_SR <- read.csv("/Users/tanishi/FYP/feature_extraction/SR/SR_binary_GO.csv")
binary_go_LR <- read.csv("/Users/tanishi/FYP/feature_extraction/LR/LR_binary_GO.csv")
binary_go_NEG <- binary_go

go_terms_SD <- colnames(binary_go_SD)[grepl("^GO_", colnames(binary_go_SD))]
go_terms_LD <- colnames(binary_go_LD)[grepl("^GO_", colnames(binary_go_LD))]
go_terms_LR <- colnames(binary_go_LR)[grepl("^GO_", colnames(binary_go_LR))]
go_terms_SR <- colnames(binary_go_SR)[grepl("^GO_", colnames(binary_go_SR))]
go_terms_NEG <- colnames(binary_go_NEG)[grepl("^GO_", colnames(binary_go_NEG))]

goslim_terms_SD <- colnames(binary_go_SD)[grepl("^GOSlim_", colnames(binary_go_SD))]
goslim_terms_SR <- colnames(binary_go_SR)[grepl("^GOSlim_", colnames(binary_go_SR))]
goslim_terms_LD <- colnames(binary_go_LD)[grepl("^GOSlim_", colnames(binary_go_LD))]
goslim_terms_LR <- colnames(binary_go_LR)[grepl("^GOSlim_", colnames(binary_go_LR))]
goslim_terms_NEG <- colnames(binary_go_NEG)[grepl("^GOSlim_", colnames(binary_go_NEG))]

# All go and goslim terms 

all_go_terms <- Reduce(union, list(go_terms_SD, go_terms_SR, go_terms_LD, go_terms_LR,
                                   go_terms_NEG))
length(all_go_terms) # n = 9233

all_go_slim_terms <- Reduce(union, list(
  goslim_terms_LD, goslim_terms_LR, goslim_terms_SD, goslim_terms_SR, goslim_terms_NEG
))

length(all_go_slim_terms) # n = 142

# Convert binary GO matrices to focus on GOSlim terms only. 

# Add missing columns (GO terms) to each set and fill with 0
for (term in setdiff(all_go_terms, colnames(binary_go_SD))) {
  binary_go_SD[[term]] <- 0
}
for (term in setdiff(all_go_terms, colnames(binary_go_LD))) {
  binary_go_LD[[term]] <- 0
}

# Reorder columns to match
binary_go_SD <- binary_go_SD[, all_go_terms]
binary_go_LD <- binary_go_LD[, all_go_terms]
# Add missing columns (GO terms) to each set and fill with 0
for (term in setdiff(all_go_terms, colnames(binary_go_SD))) {
  binary_go_SD[[term]] <- 0
}
for (term in setdiff(all_go_terms, colnames(binary_go_LD))) {
  binary_go_LD[[term]] <- 0
}

# Reorder columns to match
binary_go_SD <- binary_go_SD[, all_go_terms]
binary_go_LD <- binary_go_LD[, all_go_terms]


