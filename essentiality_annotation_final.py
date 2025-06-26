# ESSENTIALITY ANNOTATION
# --------------------------------------

# -- Load packages
import pandas as pd
import numpy as np
import itertools

# --- Load all datasets 
genotype_phenotype = pd.read_csv('/Users/tanishi/FYP/data/genotype_phenotype_data_fb_2025_02.csv') # 394694 rows x 7 columns
allele_map = pd.read_csv('/Users/tanishi/FYP/data/fbal_to_fbgn_fb_2025_02.csv') # 302333 rows x 4 columns 
unique_isoforms = pd.read_csv('/Users/tanishi/FYP/data/dmel_unique_protein_isoforms_fb_2025_02.csv') # 22459 rows x 4 columns 

# --- Explore dataset structure and organisation 
unique_isoforms.head()
genotype_phenotype.head()
allele_map.head()

# --- Ensure consistent column names 
unique_isoforms.columns = ['FBgn_ID', 'GeneSymbol', 'representative_protein', 'identifical_protein(s)']
genotype_phenotype.columns = ['AlleleSymbol', 'FBal_ID', 'phenotype_name', 'phenotype_id', 'qualifier_names', 'qualifier_ids', 'reference']
allele_map.columns = ['FBal_ID', 'AlleleSymbol', 'FBgn_ID', 'GeneSymbol']

# --- Step 1: Extract protein-coding genes with allele data available 
# Process unique isoforms and extract protein-coding gene IDs
unique_isoforms = unique_isoforms.dropna(subset = ['FBgn_ID', 'GeneSymbol']) # Remove rows with no gene ID
unique_isoforms['FBgn_ID'].nunique() # 13985 D.mel protein-coding genes in the recent FB release 
protein_coding_gene_ids = set(unique_isoforms['FBgn_ID'])

# Extract allele IDs for protein-coding genes only using unique_isoforms and allele_map 
coding_allele_map = allele_map[allele_map['FBgn_ID'].isin(protein_coding_gene_ids)]
coding_allele_map['FBgn_ID'].nunique() # n = 13500 protein-coding genes with corresponding allele IDs 
coding_allele_map.head()

# --- Step 2: Filter for alleles with phenotype annotations available 
# Clean up genotype_phenotype data by dropping alleles with missing phenotype
genotype_phenotype = genotype_phenotype.dropna(subset=['FBal_ID', 'phenotype_name'])
genotype_phenotype['FBal_ID'].nunique() # 159345 alleles with phenotypes annotated

# Remove ambiguous terms (semi-lethal, conditional, partial)
ambiguous_terms = 'semi|partial|sensitive|conditional'
mask_clear = ~genotype_phenotype['phenotype_name'].str.contains(ambiguous_terms, case=False, na=False)
mask_qualifier = ~genotype_phenotype['qualifier_names'].str.contains(ambiguous_terms, case=False, na=False)
genotype_phenotype_filtered = genotype_phenotype[mask_clear & mask_qualifier]
genotype_phenotype_filtered['FBal_ID'].nunique() # 156265 alleles with unambiguous phenotypes annotated. 3080 alleles dropped from previous step. 

# --- Step 3:  Merge phenotype annotations with protein-coding gene map
# Merge by protein-coding allele and gene IDs 
allele_phenotype_gene = genotype_phenotype_filtered.merge(
    coding_allele_map[['FBal_ID', 'FBgn_ID']],
    on='FBal_ID', # merging by Allele ID  because both data frames have this information
    how='inner'
)

# Reorder to make FBgn_ID the first column
cols = allele_phenotype_gene.columns.tolist()
cols.insert(0, cols.pop(cols.index('FBgn_ID')))
allele_phenotype_gene = allele_phenotype_gene[cols]
allele_phenotype_gene.head()

# Check structure 
print(allele_phenotype_gene.shape) # 
print(allele_phenotype_gene['FBgn_ID'].nunique()) # 5845 protein-coding genes with mapped phenotype
print(allele_phenotype_gene['FBal_ID'].nunique()) # 37576 unique protein-coding alleles with annotated phenotypes
print(allele_phenotype_gene.head())

# --- Step 4: Extract essential genes 
# List of essential phenotype keywords
essential_phenotypes = ['lethal', 'sterile']

# List of phenotype qualifiers to consider
essential_qualifiers = ['dominant', 'recessive']

# Filter alleles that have essential phenotypes
essential_rows = allele_phenotype_gene[
    allele_phenotype_gene['phenotype_name'].str.contains('|'.join(essential_phenotypes), case=False, na=False) &
    allele_phenotype_gene['qualifier_names'].str.contains('|'.join(essential_qualifiers), case=False, na=False)
]

essential_rows['FBgn_ID'].nunique() # 2081 essential genes in total, which have a lethal and/or sterile phenotype and a labelled qualifier name 
essential_rows['FBal_ID'].nunique() # 13027 essential gene alleles with labelled phenotypes. 

# Get unique essential gene IDs
essential_genes = set(essential_rows['FBgn_ID'].unique())

# --- Step 5: Extract non-essential genes. 
# Get all protein-coding genes with corresponding allele and phenotype annotations
all_protein_coding_genes = set(allele_phenotype_gene['FBgn_ID'])
non_essential_genes = all_protein_coding_genes - essential_genes
len(non_essential_genes) # 3764 non-essential genes 

# Extract all rows where the gene is non-essential
non_essential_rows = allele_phenotype_gene[allele_phenotype_gene['FBgn_ID'].isin(non_essential_genes)]

# Check number of unique non-essential genes and alleles
print(non_essential_rows['FBgn_ID'].nunique())  # 3764 protein-coding genes with non-lethal and non-sterile phenotypes
print(non_essential_rows['FBal_ID'].nunique())  # 10604 alleles 
non_essential_rows.head()

# --- Step 6: Tracking gene counts:
# all protein-coding genes with corresponding allele labels
all_protein_coding_genes = set(coding_allele_map['FBgn_ID']) 
len(all_protein_coding_genes) # 13500 genes

# all protein-coding genes with allele-phenotype annotation (from allele_phenotype_gene)
genes_with_phenotype = set(allele_phenotype_gene['FBgn_ID'])
len(genes_with_phenotype) # 5845 genes 

# all protein-coding genes classified as essential
essential_genes = set(essential_rows['FBgn_ID'].unique())
len(essential_genes) # 2081 genes

# all protein-coding genes classified as non-essential
non_essential_genes = genes_with_phenotype - essential_genes
len(non_essential_genes) # 3764 genes

# Unclassified: protein-coding genes with NO phenotype annotations
unclassified_genes = all_protein_coding_genes - genes_with_phenotype
len(unclassified_genes) # 7655 unclassified genes 

 
# --- Step 7: Sub-classify essential genes and track gene counts 
# Lethal dominant (n = 36)
lethal_dom = essential_rows[
    essential_rows['phenotype_name'].str.contains('lethal', case=False, na=False) &
    essential_rows['qualifier_names'].str.contains('dominant', case=False, na=False)]

lethal_dom['FBgn_ID'].nunique() # 36 protein-coding lethal dom genes 

# Lethal recessive (n = 1923) 
lethal_rec = essential_rows[
    essential_rows['phenotype_name'].str.contains('lethal', case=False, na=False) &
    essential_rows['qualifier_names'].str.contains('recessive', case=False, na=False)]

lethal_rec['FBgn_ID'].nunique() # 1923 protein-coding lethal dom genes 

# Sterile dominant (n = 28)
sterile_dom = essential_rows[
    essential_rows['phenotype_name'].str.contains('sterile', case=False, na=False) &
    essential_rows['qualifier_names'].str.contains('dominant', case=False, na=False)]

sterile_dom['FBgn_ID'].nunique() # 38 protein-coding lethal dom genes 

# Sterile recessive (n = 400) 
sterile_rec = essential_rows[
    essential_rows['phenotype_name'].str.contains('sterile', case=False, na=False) &
    essential_rows['qualifier_names'].str.contains('recessive', case=False, na=False)]

sterile_rec['FBgn_ID'].nunique() # 1923 protein-coding lethal dom genes 

# Accounting for gene pleiotropy -- a single 'essential' gene can have both sterile and lethal phenotypes. 
# to account for this biological complexity, I should attempt to build a multi-label classifier (OnevsRest classifier, can implement LR and SVM)
subclass_genes = set(lethal_dom['FBgn_ID']) | set(lethal_rec['FBgn_ID']) | \
                 set(sterile_dom['FBgn_ID']) | set(sterile_rec['FBgn_ID'])

len(subclass_genes)  # Should be ≤ len(essential_genes)

# Multi-label matrix 
multilabel_df = pd.DataFrame(index=list(essential_genes))
multilabel_df['lethal'] = multilabel_df.index.isin(lethal_dom['FBgn_ID']) | multilabel_df.index.isin(lethal_rec['FBgn_ID'])
multilabel_df['sterile'] = multilabel_df.index.isin(sterile_dom['FBgn_ID']) | multilabel_df.index.isin(sterile_rec['FBgn_ID'])
multilabel_df['dominant'] = multilabel_df.index.isin(lethal_dom['FBgn_ID']) | multilabel_df.index.isin(sterile_dom['FBgn_ID'])
multilabel_df['recessive'] = multilabel_df.index.isin(lethal_rec['FBgn_ID']) | multilabel_df.index.isin(sterile_rec['FBgn_ID'])

multilabel_df['essential'] = 1

# Convert booleans to integers
multilabel_df = multilabel_df.astype(int)

# Create multi-label matrix for non-essential genes: all zeros for phenotype columns
multilabel_nonessential = pd.DataFrame(index=list(non_essential_genes))
for col in ['lethal', 'sterile', 'dominant', 'recessive', 'essential']:
    multilabel_nonessential[col] = 0

# Combine both matrices
final_multilabel_df = pd.concat([multilabel_df, multilabel_nonessential])

# Optional: check shape and counts
print(final_multilabel_df.shape) # 5845 for all genes
print(final_multilabel_df['essential'].value_counts()) # 3764 non-essential, 2081 essential
print(final_multilabel_df.sum())  # sums per phenotype label. 

# --- Step 8: Save all essential vs non-essential gene lists and IDs for label matrix! 
# Convert sets to sorted lists (optional but neat)
essential_gene_list = sorted(list(essential_genes))
nonessential_gene_list = sorted(list(non_essential_genes))

# Save to desktop (update the path if needed)
essential_path = '/Users/tanishi/FYP/FlyBase_IDs/essential_genes.txt'  # change to your local desktop path if running locally
nonessential_path = '/Users/tanishi/FYP/FlyBase_IDs/nonessential_genes.txt' 

# Save gene IDs as plain text files (one gene per line)
with open(essential_path, 'w') as f:
    for gene_id in essential_gene_list:
        f.write(gene_id + '\n')

with open(nonessential_path, 'w') as f:
    for gene_id in nonessential_gene_list:
        f.write(gene_id + '\n')

print(f"Essential genes saved: {len(essential_gene_list)} to {essential_path}")
print(f"Non-essential genes saved: {len(nonessential_gene_list)} to {nonessential_path}")
