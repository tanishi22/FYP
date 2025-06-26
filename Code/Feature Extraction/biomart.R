# Extracting gene-level features with biomaRt 

install.packages("biomaRt")
library(biomaRt)
library(dplyr)
library(tidyr)
library(data.table)

# Obtain Ensembl Genomes for D. mel, linked to FlyBase DB 
ensembl <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")

#If ensembl server is down, try alternative:
ensembl <- useEnsembl(
  biomart = "ensembl",
  dataset = "dmelanogaster_gene_ensembl",
  mirror = "asia"
)

# A total of 2642 attributes can be queried from D.mel Ensembl data. 
# feature_page, structure, homology, and sequences can be queried
total_attributes <- listAttributes(ensembl) 

# Filter to obtain attributes lists for exploration
feature_attributes <- total_attributes %>%
  filter(page == "feature_page") #142 features

structure_attributes <- total_attributes %>%
  filter(page == "structure") #34 attributes

sequence_attributes <- total_attributes %>%
  filter(page == "sequences") #54 attributes

homology_attributes <- total_attributes %>%
  filter(page == "homologs") #2412 attributes

# Check filters. use flybase_gene_id to query Ensembl data for extracted FlyBase genes 
filters <- listFilters(ensembl)

# Load sterile_dom FlyBase IDs of protein-coding genes (n = 29)

# Query biomart for structure attributes (dimensions = 634 x 34)
structure_result <- getBM(
  attributes = structure_attributes$name,
  filters = "flybase_gene_id",
  values = NEG_4_gene_isoform$FBgn_ID, # Refer to gene_isoform_mapping.R to access this! 
  mart = ensembl
) 

# Note, this returns transcript- and protein isoform-specific features for each unique gene ID
# I will extract the structural features for the representative protein isoforms for each gene, 
# which was extracted with a separate Python script (see isoform_processing.py)

# Filter structure_result to retain only rows with matching FBpp IDs (dim = 150 x 34)
filtered_structure <- structure_result[structure_result$ensembl_peptide_id 
                                       %in% NEG_4_gene_isoform$ensembl_peptide_id, ]
filtered_structure_copy <- filtered_structure

# For aggregation, average across numeric features and group for categorical features

# First, isolate numeric columns only
numeric_cols <- filtered_structure %>%
  select(where(is.numeric)) %>%
  colnames()

# Then group and summarise across those
numeric_structure_agg <- filtered_structure %>%
  group_by(ensembl_peptide_id) %>%
  summarise(across(all_of(numeric_cols), mean, na.rm = TRUE), .groups = "drop")

# add protein isoform count length to numeric dataset
isoform_count <- NEG_4_isoforms[, c("Num_Isoforms", 'ensembl_peptide_id', 'Longest_Isoform_Length')]
numeric_structure_agg <- left_join(numeric_structure_agg, isoform_count, by = 'ensembl_peptide_id')

# aggregated categorical features by mode rather than grouping so that the data is more compatible 
# for downstream one-hot encoding and machine learning. but unclear if this is the best approach? (Apr 28)
# after filtering on isoform-level, this makes no difference so it's fine! (May 6)
categorical_structure_agg <- filtered_structure %>%
  group_by(ensembl_peptide_id) %>%
  summarise(across(where(is.character), 
                   ~names(sort(table(.), decreasing = TRUE))[1], 
                   .names = "mode_{.col}")) 


# Concatenate numerical and categorical structure features (29 x 34 data table)
structure_features <- left_join(numeric_structure_agg, categorical_structure_agg, by = "ensembl_peptide_id")

# Remove irrelevant features (eg features that show no variance between gene labels)
# ie external_gene_source is 'gene_name' for all genes. 
# transcript_start and _end capture the same info as start_ and end_position so removed
structure_features <- structure_features[, !names(structure_features) 
                                         %in% c('mode_external_gene_source', 
                                                'mode_gene_biotype',
                                                'transcript_start',
                                                'transcript_end',
                                                'mode_ensembl_transcript_id',
                                                'mode_external_gene_name',
                                                'mode_description',
                                                'mode_ensembl_exon_id')] 

# structure_features has gene, transcript, and peptide IDs 
# overall, 28 structural features extracted for 29 genes 

# add features like gene length, transcript length, etc 
setDT(structure_features)
structure_features[, gene_length := abs(end_position - start_position) + 1]
structure_features[, coding_gene_length := abs(genomic_coding_end - genomic_coding_start) + 1]
structure_features[, five_utr_length := abs(`5_utr_end` - `5_utr_start`) + 1]
structure_features[, three_utr_length := abs(`3_utr_end` - `3_utr_start`) + 1]
structure_features[, cdna_coding_length := abs(cdna_coding_end - cdna_coding_start) + 1]
structure_features[, exon_chrom_length := abs(exon_chrom_end - exon_chrom_start) + 1]

summary(structure_features) # 32 features
write.csv(structure_features, "~/FYP/feature_extraction/structure/NEG_4_structure_features.csv", 
          row.names = FALSE)
# Repeat for domain features, gene ontology features, and GC content from feature_attributes

# Domain features 
domain_attributes <- c('ensembl_gene_id','hmmpanther', 'pfam', 'prints', 'scanprosite',
                       'pfscan', 'smart', 'superfamily')

# Secondary structure features
ss_attributes <- c('ensembl_gene_id','ncoils', 'seg','signalp','tmhmm')

# Extract domain features (764 x 8)
domain_result <- getBM(
  attributes = c('ensembl_gene_id', 'pfam', 'prints', 'scanprosite', 'pfscan', 'smart', 
                 'superfamily'),
  filters = "flybase_gene_id",
  values = EG_1_gene_isoform$FBgn_ID,
  mart = ensembl
) 

# Extract secondary structure features (34 x 5) majority of this is missing, so will need to 
# use other predictie methods to get this data, as done in other papers
ss_result <- getBM(
  attributes = ss_attributes,
  filters = "flybase_gene_id",
  values = lethal_dom_ids,
  mart = ensembl
) 

# Domain results can be processed readily, as we're interested in if a gene as a domain attribute annotaiton
# not WHAT annotation it has. The latter is to specific and sparse, which will unnecessarily 
# increase the feature matrix size. 

domain_features <- domain_result %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    # Check if there's *any non-empty* annotation in each domain column
    has_pfam = any(!is.na(pfam) & pfam != ""),
    has_hmmpanther = any(!is.na(hmmpanther) & hmmpanther != ""),
    has_smart = any(!is.na(smart) & smart != ""),
    has_scanprosite = any(!is.na(scanprosite) & scanprosite != ""),
    has_superfamily = any(!is.na(superfamily) & superfamily != ""),
    has_pfscan = any(!is.na(pfscan) & pfscan != ""),
    has_prints = any(!is.na(prints) & prints != ""),
    .groups = "drop"
  ) %>%
  # Convert logical TRUE/FALSE to binary 1/0
  mutate(across(starts_with("has_"), ~ as.integer(.)))

write.csv(domain_features, "~/FYP/feature_extraction/SR/SR_domain_features.csv", 
          row.names = TRUE)
