# Load package and data. Starting with sterile dom protein sequences. 
library(protr)

cds <- readFASTA('~/FYP/prot_sequences/NEG_translation.fasta')
length(cds) # n = 29

# Check if protein sequences have standard amino acids and process. 
cds <- cds[(sapply(cds, protcheck))] 
length(cds) # n = 28. one sequence removed as it did not have standard amino acids, revise which one was excluded! 

# Amino acid composition features (n = 420). 
# Note, can choose to exclude tripeptide comp as this descriptor has 8000 features
aac <- sapply(cds, protr::extractAAC) # AA comp (n=20)
t_aac <- t(aac)
dc <- t(sapply(cds, protr::extractDC)) # Dipeptide comp (n = 400)
tc <- t(sapply(cds, protr::extractTC)) # Tripeptide comp (n = 8000)
aa_comp <- data.frame(t_aac, dc)
aa_tc_comp <- data.frame(aa_comp, tc) # n = 8420
write.csv(aa_tc_comp, "~/FYP/feature_extraction/NEG/NEG_aa_tc_comp.csv", row.names = FALSE)

# For NEGs, need to split this up into batches!!
# Helper function to split into batches
split_batches <- function(x, batch_size) {
  split(x, ceiling(seq_along(x) / batch_size))
}

# Set batch size (adjust based on memory/performance)
batch_size <- 500
cds_batches <- split_batches(cds, batch_size)

# Apply feature extraction in batches
extract_feature_batch <- function(batches, extractor_func, ...) {
  do.call(rbind, lapply(batches, function(batch) {
    t(sapply(batch, extractor_func, ...))
  }))
}


# Extract autocorrelation features (n = 720), 240 features for each autocorrelation method 
?extractMoreauBroto
moreau <- t(sapply(cds, protr::extractMoreauBroto)) 
moran <- t(sapply(cds, protr::extractMoran)) 
geary <- t(sapply(cds, protr::extractGeary))
autocorrelation <- data.frame(moreau, moran, geary) 

# Autocorrelation features in batches
moreau <- extract_feature_batch(cds_batches, protr::extractMoreauBroto)
moran  <- extract_feature_batch(cds_batches, protr::extractMoran)
geary  <- extract_feature_batch(cds_batches, protr::extractGeary)
autocorrelation <- data.frame(moreau, moran, geary)

# Extract CTD features (n = 147) - Global distribution of AA properties, summarising physicochemical trends over AAseqs
ctdc <- t(sapply(cds, protr::extractCTDC)) # Composition, n = 21 
ctdt <- t(sapply(cds, protr::extractCTDT)) # Transition, n = 21 
ctdd <- t(sapply(cds, protr::extractCTDD)) # Distribution, n = 105
ctd <- data.frame(ctdc, ctdt, ctdd)

ctdc <- extract_feature_batch(cds_batches, protr::extractCTDC)
ctdt <- extract_feature_batch(cds_batches, protr::extractCTDT)
ctdd <- extract_feature_batch(cds_batches, protr::extractCTDD)
ctd  <- data.frame(ctdc, ctdt, ctdd)

# Conjoint triad features (n = 343)
ctriad <- t(sapply(cds, protr::extractCTriad)) 
ctriad <- extract_feature_batch(cds_batches, protr::extractCTriad)

# Quasi-sequence order features (n =  100)
qso <- t(sapply(cds, protr::extractQSO, nlag = 10)) # Quasi-sequence-order descriptor
min(sapply(cds, nchar))
qso <- extract_feature_batch(cds_batches, protr::extractQSO, nlag = 10)

# Pseudo amino acid composition (n = 50)
PseAAC <- t(sapply(cds, protr::extractPAAC)) # Pseudo-AA comp (n = 50)
APseAAC <- t(sapply(cds, protr::extractAPAAC)) # Amphiphilic pseudo-AA comp (n = 80)
pseudo_AA <- data.frame(PseAAC, APseAAC)
??PseAAC
PseAAC  <- extract_feature_batch(cds_batches, protr::extractPAAC)
APseAAC <- extract_feature_batch(cds_batches, protr::extractAPAAC, lambda = 10)
pseudo_AA <- data.frame(PseAAC, APseAAC)

#merge all matrices into single feature matrix (28 x 1860)
NEG_protr <-  data.frame(autocorrelation, ctd, ctriad, qso, pseudo_AA)

write.csv(autocorrelation, "~/FYP/feature_extraction/NEG/NEG_autocorrelation.csv", 
          row.names = FALSE)

write.csv(ctriad, "~/FYP/feature_extraction/NEG/NEG_ctd_protr.csv", 
          row.names = FALSE)

# add peptide IDs and gene ids
NEG_protr$ensembl_peptide_id <- rownames(NEG_protr)
rownames(NEG_protr) <- NULL
NEG_protr_final <- left_join(NEG_gene_isoform, NEG_protr, by = 'ensembl_peptide_id')

# save
write.csv(NEG_protr_final, "~/FYP/feature_extraction/NEG/NEG_protr_matrix.csv", row.names = FALSE)

