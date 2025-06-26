# Install packages
library(seqinr)
library(protr)

# Load files with seqinr 
cds_seqin <- read.fasta("~/FYP/CDS/processed/sterile_dom_isoform.fasta", 
                        seqtype = 'DNA', as.string = FALSE) #nt sequences 
prot_seqin <- read.fasta('~/FYP/prot_sequences/SD_translation.fasta', 
                         seqtype = "AA", as.string = FALSE) #FB-curated aa sequences for representative isoform

gene_isoform_map <- read.csv("~/FYP/feature_extraction/SD/sterile_dom_isoforms.csv")
gene_isoform_map <- gene_isoform_map[, c("FBgn_ID", "Longest_Isoform_ID")]
colnames(gene_isoform_map)[colnames(gene_isoform_map) == 'Longest_Isoform_ID'] <- 'ensembl_peptide_id'

# Starting with seqinr feature extraction:
# aa features of interest: pI, physicochemical classes, aa counts, codon adaptation index 
# nt features of interest: GC content

gc_content <- sapply(cds_seqin, GC) #extract GC content for each isoform
theoretical_pi <- sapply(prot_seqin, seqinr::computePI)

seqinr_features <- data.frame(gc_content, theoretical_pi)
seqinr_features$ensembl_peptide_id <- rownames(seqinr_features)
seqinr_features <- cbind(ensembl_peptide_id, seqinr_features)
rownames(seqinr_features) <- NULL 

# Extract proportion of physicochemical classes
# Amino acid composition extracted from protr and theoretical pI extracted from seqinr computePI()
# so focus on physicochem class proportions 
aa_stats <- sapply(prot_seqin, AAstat)
prop_list <- lapply(aa_stats["Prop", ], function(x) as.data.frame(t(unlist(x))))
prop_df <- do.call(rbind, prop_list)
prop_df$ensembl_peptide_id <- rownames(prop_df)
rownames(prop_df) <- NULL

# merge dataframes 
seqinr_features <- left_join(seqinr_features, prop_df, by = 'ensembl_peptide_id')
seqinr_features <- left_join(gene_isoform_map, seqinr_features, by = 'ensembl_peptide_id')

write.csv(seqinr_features, "~/FYP/feature_extraction/SD/SD_seqinr.csv")
