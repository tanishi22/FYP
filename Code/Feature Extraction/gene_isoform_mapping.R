library(tidyr)
library(tidyverse)
library(dplyr)

# extract isoform maps for egs and negs 
getwd()

# EG 1
EG_1_isoforms <- read.csv("~/FYP/CDS/processed/EG_1_isoforms.csv")
EG_1_isoforms$X = NULL
colnames(EG_1_isoforms)[colnames(EG_1_isoforms) == "Longest_Isoform_ID"] <- "ensembl_peptide_id"
EG_1_gene_isoform <- EG_1_isoforms[, c("FBgn_ID", "ensembl_peptide_id")]

# EG 2
EG_2_isoforms <- read.csv("~/FYP/CDS/processed/EG_2_isoforms.csv")
EG_2_isoforms$X = NULL
colnames(EG_2_isoforms)[colnames(EG_2_isoforms) == "Longest_Isoform_ID"] <- "ensembl_peptide_id"
EG_2_gene_isoform <- EG_2_isoforms[, c("FBgn_ID", "ensembl_peptide_id")]

# NEG 1
NEG_1_isoforms <- read.csv("~/FYP/CDS/processed/NEG_1_isoforms.csv")
NEG_1_isoforms$X = NULL
colnames(NEG_1_isoforms)[colnames(NEG_1_isoforms) == "Longest_Isoform_ID"] <- "ensembl_peptide_id"
NEG_1_gene_isoform <- NEG_1_isoforms[, c("FBgn_ID", "ensembl_peptide_id")]

# NEG 2
NEG_2_isoforms <- read.csv("~/FYP/CDS/processed/NEG_2_isoforms.csv")
NEG_2_isoforms$X = NULL
colnames(NEG_2_isoforms)[colnames(NEG_2_isoforms) == "Longest_Isoform_ID"] <- "ensembl_peptide_id"
NEG_2_gene_isoform <- NEG_2_isoforms[, c("FBgn_ID", "ensembl_peptide_id")]

# NEG 3
NEG_3_isoforms <- read.csv("~/FYP/CDS/processed/NEG_3_isoforms.csv")
NEG_3_isoforms$X = NULL
colnames(NEG_3_isoforms)[colnames(NEG_3_isoforms) == "Longest_Isoform_ID"] <- "ensembl_peptide_id"
NEG_3_gene_isoform <- NEG_3_isoforms[, c("FBgn_ID", "ensembl_peptide_id")]

# NEG 4
NEG_4_isoforms <- read.csv("~/FYP/CDS/processed/NEG_4_isoforms.csv")
NEG_4_isoforms$X = NULL
colnames(NEG_4_isoforms)[colnames(NEG_4_isoforms) == "Longest_Isoform_ID"] <- "ensembl_peptide_id"
NEG_4_gene_isoform <- NEG_4_isoforms[, c("FBgn_ID", "ensembl_peptide_id")]

# Lethal dom
LD_gene_isoform <- read.csv('~/FYP/feature_extraction/LD/lethal_dom_isoforms.csv')
LD_gene_isoform$X <- NULL
colnames(LD_gene_isoform)[colnames(LD_gene_isoform) == "Longest_Isoform_ID"] <- "ensembl_peptide_id"
LD_gene_isoform <- LD_gene_isoform[, c('FBgn_ID', 'ensembl_peptide_id' )]

lethal_dom_ids <- LD_gene_isoform$FBgn_ID
LD_longest_ids <- LD_gene_isoform$ensembl_peptide_id


# Lethal rec
LR_gene_isoform <- read.csv('~/FYP/feature_extraction/LR/lethal_rec_isoforms.csv')
LR_gene_isoform$X <- NULL
colnames(LR_gene_isoform)[colnames(LR_gene_isoform) == "Longest_Isoform_ID"] <- "ensembl_peptide_id"
LR_gene_isoform <- LR_gene_isoform[, c('FBgn_ID', 'ensembl_peptide_id' )]

lethal_rec_ids <- LR_gene_isoform$FBgn_ID
LR_longest_ids <- LR_gene_isoform$ensembl_peptide_id

# Sterile dom
SD_gene_isoform <- read.csv('~/FYP/feature_extraction/SD/sterile_dom_isoforms.csv')
SD_gene_isoform$X <- NULL
colnames(SD_gene_isoform)[colnames(SD_gene_isoform) == "Longest_Isoform_ID"] <- "ensembl_peptide_id"
SD_gene_isoform <- SD_gene_isoform[, c('FBgn_ID', 'ensembl_peptide_id' )]

sterile_dom_ids <- SD_gene_isoform$FBgn_ID
SD_longest_ids <- SD_gene_isoform$ensembl_peptide_id

# Sterile rec
SR_gene_isoform <- read.csv('~/FYP/feature_extraction/SR/sterile_rec_isoforms.csv')
SR_gene_isoform$X <- NULL
colnames(SR_gene_isoform)[colnames(SR_gene_isoform) == "Longest_Isoform_ID"] <- "ensembl_peptide_id"
SR_gene_isoform <- SR_gene_isoform[, c('FBgn_ID', 'ensembl_peptide_id' )]

sterile_rec_ids <- SR_gene_isoform$FBgn_ID
SR_longest_ids <- SR_gene_isoform$ensembl_peptide_id

# Non-essential genes 
NEG_gene_isoform <- read.csv('~/FYP/feature_extraction/NEG/neg_isoforms.csv')
NEG_gene_isoform$X <- NULL
colnames(NEG_gene_isoform)[colnames(NEG_gene_isoform) == "Longest_Isoform_ID"] <- "ensembl_peptide_id"
NEG_gene_isoform <- NEG_gene_isoform[, c('FBgn_ID', 'ensembl_peptide_id' )]

NEG_ids <- NEG_gene_isoform$FBgn_ID
NEG_longest_ids <- NEG_gene_isoform$ensembl_peptide_id