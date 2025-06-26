# Libraries
library(data.table)
library(tidyr)
library(tidyverse)

# Load data 
paralogs_data <- fread("~/FYP/data/homology/dmel_paralogs_fb_2025_02.tsv", skip = 4)

# Start with paralogs analysis (223056 obs x 11 variables)
# Features of interest: 
# Number of paralogs per gene, DIOPT score summary stats (mean, median, max, min), 
# distance from gene if on same chromosomal arm, 
# Categorical: presence/absence of paralogs, presence/absence on same arm/scaffold

# First, organise raw dataset 
colnames(paralogs_data)[colnames(paralogs_data) == "## FBgn_ID"] <- "FBgn_ID"
colnames(paralogs_data)[colnames(paralogs_data) == "Arm/Scaffold"] <- "Arm_Scaffold"
colnames(paralogs_data)[colnames(paralogs_data) == "Paralog_Arm/Scaffold"] <- "Paralog_Arm_Scaffold"
paralogs_data <- paralogs_data[, .(FBgn_ID, Arm_Scaffold, Location, Paralog_FBgn_ID,
                                   Paralog_Arm_Scaffold, Paralog_Location, DIOPT_score)]

# load gene IDs
EG_ids <- read.table('~/FYP/Flybase_IDs/essential_genes.txt')
NEG_ids <- read.table('~/FYP/Flybase_IDs/nonessential_genes.txt')
# Subset based on essential gene IDs (refer to gene set-specific isoform map)
essential_paralogs <- paralogs_data[FBgn_ID %in% EG_ids$V1] # 892 x 7
nonessential_paralogs <- paralogs_data[FBgn_ID %in% NEG_ids$V1]

# Calculate number of unique paralogs per essential genes. 
paralog_counts <- nonessential_paralogs %>%
  group_by(FBgn_ID) %>%
  summarise(Num_Paralogs = n_distinct(Paralog_FBgn_ID))

# Note, only 25 genes have paralogs so will need to add 0 for the rest 
missing_genes <- setdiff(NEG_ids$V1, paralog_counts$FBgn_ID) # n = 4

# Create a df for these missing genes and set their paralog count to 0
missing_paralog_counts <- data.frame(FBgn_ID = missing_genes, Num_Paralogs = 0)

# Combine the missing genes with the existing paralog counts and sort by FBgn_ID for consistency
paralog_counts <- bind_rows(paralog_counts, missing_paralog_counts)
paralog_counts <- paralog_counts %>% arrange(FBgn_ID)

# Presence/absence of paralogs
paralog_counts <- paralog_counts %>%
  mutate(paralogs_presence = ifelse(Num_Paralogs > 0, 1, 0))

# Calculate mean, median, max, and min DIOPT score for each gene
diopt_summary <- nonessential_paralogs %>%
  group_by(FBgn_ID) %>%
  summarise(
    para_mean_diopt_score = mean(DIOPT_score, na.rm = TRUE),
    para_median_diopt_score = median(DIOPT_score, na.rm = TRUE),
    para_max_diopt_score = max(DIOPT_score, na.rm = TRUE),
    para_min_diopt_score = min(DIOPT_score, na.rm = TRUE)
  )

# Merge the DIOPT summary with the paralog_counts data
paralog_counts <- paralog_counts %>%
  left_join(diopt_summary, by = "FBgn_ID")

# For genes without any paralogs (i.e., where paralogs_presence is 0), set DIOPT scores to 0 or NA
paralog_counts <- paralog_counts %>%
  mutate(
    para_mean_diopt_score = ifelse(paralogs_presence == 0, 0, para_mean_diopt_score),
    para_median_diopt_score = ifelse(paralogs_presence == 0, 0, para_median_diopt_score),
    para_max_diopt_score = ifelse(paralogs_presence == 0, 0, para_max_diopt_score),
    para_min_diopt_score = ifelse(paralogs_presence == 0, 0, para_min_diopt_score)
  )

# Create a new variable indicating if any paralogs are on the same arm/scaffold
paralog_same_arm <- nonessential_paralogs %>%
  group_by(FBgn_ID) %>%
  summarise(
    paralog_same_arm = as.integer(any(Arm_Scaffold == Paralog_Arm_Scaffold))
  )

# Merge with paralog_counts or your existing data
paralog_counts <- paralog_counts %>%
  left_join(paralog_same_arm, by = "FBgn_ID")


# Calculate the average distance between each gene and it's paralog
# First, extract start and end locations and calculate the distance 
# Convert location format to numeric start and end
nonessential_paralogs <- nonessential_paralogs %>%
  mutate(
    Gene_Start = as.numeric(sub("^(\\d+)\\..*", "\\1", Location)),
    Gene_End = as.numeric(sub("^\\d+..(\\d+)", "\\1", Location)),
    Paralog_Start = as.numeric(sub("^(\\d+)\\..*", "\\1", Paralog_Location)),
    Paralog_End = as.numeric(sub("^\\d+..(\\d+)", "\\1", Paralog_Location))
  )

# Calculate distance only for paralogs on the same arm/scaffold as the gene
nonessential_paralogs_same_arm <- nonessential_paralogs %>%
  filter(Arm_Scaffold == Paralog_Arm_Scaffold) %>%
  mutate(
    Gene_Midpoint = (Gene_Start + Gene_End) / 2,
    Paralog_Midpoint = (Paralog_Start + Paralog_End) / 2,
    Distance = abs(Gene_Midpoint - Paralog_Midpoint)
  )


# Calculate the average distance for each gene, only considering paralogs on the same arm/scaffold
paralog_distance_avg <- nonessential_paralogs_same_arm %>%
  group_by(FBgn_ID) %>%
  summarise(
    avg_distance_same_arm = mean(Distance, na.rm = TRUE)
  )

# Merge the average distance with paralog_counts or your final feature matrix
paralog_counts <- paralog_counts %>%
  left_join(paralog_distance_avg, by = "FBgn_ID")

paralog_features <- paralog_counts
write.csv(paralog_features, "~/FYP/feature_extraction/NEG_paralog_features.csv", 
          row.names = FALSE)