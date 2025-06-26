# Libraries
library(data.table)
library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)

# ------------------------------- Species-Specific DIOPT Stats ----------------------------------
# --- Features of interest:
# No. of databases that the prediction was derived from per gene for ortholog pre-filtering 
# Presence/absence of orthologs, Total ortholog count, species-specific ortholog count per gene 
# Avg diopt score per gene for each species, diopt summary stats per gene (min, max, mean, median)
# No. of high ranking orthologs per gene 

# --- Load data 
# Note, all of this data was downloaded from the FlyRNAi DIOPT server, with the sterile dom gene ID list as input
# Need to change input and download corresponding files for each gene set!! 
SR_diopt <- read.csv('~/FYP/feature_extraction/diopt/raw/NEG3_b2_diopt.csv') # (orthology mapped to 12 model species)

# --- Exclude irrelevant information and rename columns names 
colnames(SR_diopt)
SR_diopt <- SR_diopt %>% 
  select(-`Search.Term`, -`Fly.GeneID`, -`All.GeneID`, -`All.Symbol`, -`Ensmbl.ID..link.HPA.`,
        -`Alignment...Scores`, -Feedback, -`Gene.Details`, -`Input.Order`)

colnames(SR_diopt)[colnames(SR_diopt) == "FlyBaseID"] <- "FBgn_ID"
colnames(SR_diopt)[colnames(SR_diopt) == "Species 2"] <- "Species"
# --- To ensure ortholog counts are reliable, I will only include orthology predictions/annotations 
#     if they are supported by 2 or more databases. 
#     This can be determined from the Prediction_Derived_From column 

# Count number of sources using str_count and filter out based on criteria -> 49 rows removed
SR_diopt_filtered <- SR_diopt %>%
  mutate(num_sources = str_count(Prediction.Derived.From, ",") + 1) %>%
  filter(num_sources >= 2)

# --- Calculate ortholog count per gene - 2 genes with no orthologs based on defined criteria
# Add zeros for all missing genes' values after all the features have been extracted
ortholog_counts <- SR_diopt_filtered %>%
  filter(!is.na(All.Species.Gene.ID) & All.Species.Gene.ID != "") %>%
  group_by(FBgn_ID) %>%
  summarise(Ortholog_Count = n_distinct(All.Species.Gene.ID), .groups = "drop")

# missing_orthos <- setdiff(SR_gene_isoform$FBgn_ID, ortholog_counts$FBgn_ID)
# missing_ortho_counts <- data.frame(FBgn_ID = missing_orthos, Ortholog_Count = 0)

# --- Presence/absence of orthologs
# doesn't seem like a useful feature now, but let's see if bigger gene lists show more variance in this feature. 
ortholog_counts <- ortholog_counts %>%
  mutate(ortholog_presence = ifelse(Ortholog_Count > 0, 1, 0)) 

# --- Ortholog count per gene for each species 
ortholog_counts_species <- SR_diopt_filtered %>%
  filter(!is.na(All.Species.Gene.ID) & All.Species.Gene.ID != "") %>%
  group_by(FBgn_ID, Species.2) %>%
  summarise(Species_Ortholog_Count = n_distinct(All.Species.Gene.ID), .groups = "drop")

# Flip
ortholog_counts_wide <- ortholog_counts_species %>%
  pivot_wider(names_from = Species.2,
              values_from = Species_Ortholog_Count,
              names_glue = "{Species.2}_Ortholog_Count", # Specify each species
              values_fill = 0)  # fill missing combinations with 0

# --- Mean, median, max, and min DIOPT score for each gene and 0 if N/A
diopt_gene_stats <- SR_diopt_filtered %>%
  group_by(FBgn_ID) %>%
  summarise(
    ortho_mean_diopt   = if (all(is.na(DIOPT.Score))) 0 else mean(DIOPT.Score, na.rm = TRUE),
    ortho_median_diopt = if (all(is.na(DIOPT.Score))) 0 else median(DIOPT.Score, na.rm = TRUE),
    ortho_max_diopt    = if (all(is.na(DIOPT.Score))) 0 else max(DIOPT.Score, na.rm = TRUE),
    ortho_min_diopt    = if (all(is.na(DIOPT.Score))) 0 else min(DIOPT.Score, na.rm = TRUE),
    .groups = "drop"
  )

# --- Average DIOPT score per species 
# Group by gene and species, and calculate the mean DIOPT score
diopt_species_pair_avg <- SR_diopt_filtered %>%
  group_by(FBgn_ID, Species.2) %>%
  summarise(avg_diopt_score = mean(DIOPT.Score, na.rm = TRUE), .groups = "drop") %>%
  mutate(species_pair = paste('Fly', Species.2, 'AvgDIOPT' ,sep = "_"))

# Use dcast to reshape the data into wide format
diopt_species_pair_wide <- dcast(diopt_species_pair_avg, 
                                 FBgn_ID ~ species_pair, 
                                 value.var = "avg_diopt_score", 
                                 fill = 0)

# Merge with other diopt ortholog stats 
diopt_gene_stats <- left_join(diopt_gene_stats, diopt_species_pair_wide,
                              by = "FBgn_ID")

# --- No. of high ranking orthologs per gene 
high_rank_counts <- SR_diopt_filtered %>%
  filter(Rank == "high") %>%
  group_by(FBgn_ID) %>%
  summarise(High_Rank_Ortholog_Count = n_distinct(All.Species.Gene.ID), .groups = "drop")

# Account for missing data 
# missing_high_rank <- setdiff(SR_gene_isoform$FBgn_ID, high_rank_counts$FBgn_ID)
# missing_high_counts <- data.frame(FBgn_ID = missing_high_rank, High_Rank_Ortholog_Count = 0)

# --- Presence/absence of high-ranking orthologs 
high_rank_counts <- high_rank_counts %>%
  mutate(high_ranking_ortholog_presence = ifelse(High_Rank_Ortholog_Count > 0, 1, 0)) 

# -- Concatenate
feature_tables <- list(ortholog_counts, ortholog_counts_wide, diopt_gene_stats, high_rank_counts)

orthology_features <- reduce(feature_tables, full_join, by = "FBgn_ID") %>%
  replace(is.na(.), 0) %>%
  arrange(FBgn_ID)

# --- Account for missing genes:
#missing_genes <- setdiff(sterile_rec_ids, orthology_features$FBgn_ID)
#empty_features <- as.data.frame(matrix(0, nrow = length(missing_genes), 
                                       # ncol = ncol(orthology_features) - 1))
#colnames(empty_features) <- colnames(orthology_features)[-1]
#empty_features$FBgn_ID <- missing_genes
#empty_features <- empty_features[, colnames(orthology_features)]

# Add the missing rows to your existing table
#orthology_features <- bind_rows(orthology_features, empty_features) %>%
#  arrange(FBgn_ID)

# --- Save features
write.csv(orthology_features, "~/FYP/feature_extraction/diopt/processed/NEG3_b2_orthology.csv", 
          row.names = FALSE)

# -------------------------- Human Disease Orthologs -------------------------------------- #
# --- DIOPT-DIST results, available from FB_2025_02 Current Release
# Model organism genes are mapped to putative orthologs in the human genome, and human genes 
# are associated with disease-gene information

# --- Why do we care?
# If a Dmel gene has a high-confidence human ortholog that is disease-associated, that strongly 
# suggests it's functionally important and likely essential. 
# Would be interesting to see the predictive strength of this feature 

# --- Relevant features:
# Has_Disease_Associated_Human_Ortholog
# Num_Disease_Associations

human_disease_homologs <- fread("~/FYP/data/dmel_human_orthologs_disease_fb_2025_02.tsv", 
                                skip = 3)
colnames(human_disease_homologs)[colnames(human_disease_homologs) == "##Dmel_gene_ID"] <- "FBgn_ID"
subset_disease_homologs <- human_disease_homologs[FBgn_ID %in% NEG_4_gene_isoform$FBgn_ID]

# --- Calculate homolog count per gene - 4/29 genes with no disease-associated human ortholog
disease_homolog_counts <- subset_disease_homologs %>%
  group_by(FBgn_ID) %>%
  summarise(Disease_Homolog_Count = n_distinct(Human_gene_OMIM_ID), .groups = "drop")

# Add zeros for all missing genes' values
missing_disease_homologs <- setdiff(NEG_4_gene_isoform$FBgn_ID, disease_homolog_counts$FBgn_ID)
missing_disease_counts <- data.frame(FBgn_ID = missing_disease_homologs, Disease_Homolog_Count = 0)

# Append 'missing' genes
disease_homolog_counts <- bind_rows(disease_homolog_counts, missing_disease_counts)

# save 
write.csv(disease_homolog_counts, "~/FYP/feature_extraction/diopt/NEG4_disease_homologs.csv",
          row.names = FALSE )

# --- Presence/absence of human disease orthology
# doesn't seem like a useful feature now, but let's see if bigger gene lists show more variance in this feature. 
disease_homolog_counts <- disease_homolog_counts %>%
  group_by(FBgn_ID) %>%
  mutate(has_disease_ortholog = ifelse(Disease_Homolog_Count > 0, 1, 0))

# ----------------------------- Concatenate Features -------------------------------------- # 
# --- Concatenate all orthology features with each other (note, no. of rows/genes differ in some features)
feature_tables <- list(ortholog_counts, ortholog_counts_wide, diopt_gene_stats, high_rank_counts,
                       disease_homolog_counts)

orthology_features <- reduce(feature_tables, full_join, by = "FBgn_ID") %>%
  replace(is.na(.), 0) %>%
  arrange(FBgn_ID)

# --- Account for missing genes:
missing_genes <- setdiff(sterile_rec_ids, orthology_features$FBgn_ID)
empty_features <- as.data.frame(matrix(0, nrow = length(missing_genes), 
                                       ncol = ncol(orthology_features) - 1))
colnames(empty_features) <- colnames(orthology_features)[-1]
empty_features$FBgn_ID <- missing_genes
empty_features <- empty_features[, colnames(orthology_features)]

# Add the missing rows to your existing table
orthology_features <- bind_rows(orthology_features, empty_features) %>%
  arrange(FBgn_ID)

# --- Save features
write.csv(orthology_features, "~/FYP/feature_extraction/SR/SR_orthology_features.csv", 
          row.names = FALSE)