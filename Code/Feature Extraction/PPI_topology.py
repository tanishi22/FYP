# Load packages 

import pandas as pd
import numpy as np
import os  
import subprocess
import sys
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import community as community_louvain

subprocess.check_call([sys.executable, "-m", "pip", "install", "python-louvain"])
import community as community_louvain

# --- Load and organise data 
file_path = "~/FYP/data/physical_interactions_mitab_fb_2025_02.tsv"
mitab_data = pd.read_table(file_path) 
print(mitab_data.columns)
cols_to_exclude = ['Publication 1st Author(s)', 'Publication ID(s)', 'Expansion Method(s)', 'Xref(s) Interactor A', 'Xref(s) Interactor B', 'Host Organism(s)', 'Creation Date', 'Update Date', 'Negative', 'Feature(s) Interactor A',
                   'Feature(s) Interactor B', 'Stoichiometry Interactor A', 'Stoichiometry Interactor B', 'Identification Method(s) Participant A', 
                   'Identification Method(s) Participant B', 'Confidence Value(s)', 'Checksum Interactor A', 'Checksum Interactor B', 'Interaction Checksum',
                   'Interaction Parameters', 'Biological Role(s) Interactor A', 'Biological Role(s) Interactor B', 'Source Database(s)', 'Interaction Parameters']

mitab_filtered = mitab_data.drop(columns = cols_to_exclude)
mitab_filtered.columns

# --- Extract FlyBase gene IDs (remove 'flybase:' prefix)
mitab_filtered['#ID(s) Interactor A'] = mitab_filtered.iloc[:, 0].str.replace('flybase:', '', regex=False)
mitab_filtered['ID(s) Interactor B'] = mitab_filtered.iloc[:, 1].str.replace('flybase:', '', regex=False)

# --- Load gene IDs and filter mitab
# EG_ids = pd.read_csv('/Users/tanishi/FYP/FlyBase_IDs/essential_genes.txt', header = None) # using the protein-coding gene ID list 
# EG_ids = EG_ids[0].tolist()

# mitab_EG = mitab_filtered[
    # mitab_filtered['#ID(s) Interactor A'].isin(EG_ids) |
   #  mitab_filtered['ID(s) Interactor B'].isin(EG_ids)
# ] # 37437 rows extracted

# --- Remove self-interactions
mitab_clean = mitab_filtered[
    mitab_filtered['#ID(s) Interactor A'] != mitab_filtered['ID(s) Interactor B']
]

# --- Remove self-interactions (where A == B), NetworkX will not filter this out. 26 rows filtered out. 
#mitab_EG_unique = mitab_EG[
   # mitab_EG['#ID(s) Interactor A'] != mitab_EG['ID(s) Interactor B']
#] # 36429 rows extracted. 1008 self-interactions removed. 

# --- when constructing the PPI network (PPIn), certain key assumptions have been made 
# (1) undirected graph: visualises simple connections between nodes with no flow or direction implied. this is the standard approach for PPIn analysis
# but it's also possible to add 'weights' to the edges, depending on prior literature or the type of experimental detection method used. 
# Several references to support this assumtion, so move ahead with undirected graph. 
# (2) In the network file, there shouLR be a 'source' and 'target' column. We're assuming that interactor A (column 1) is the source and interactor B (column 2) 
# is the target. 

# --- Construct graph
G = nx.from_pandas_edgelist(
    mitab_clean,
    source='#ID(s) Interactor A',
    target='ID(s) Interactor B'
)

# --- Basic overall visualisation as preliminary checkpoint  - works fine!
plt.figure(figsize=(10, 10))
nx.draw(
    G, 
    with_labels=False, 
    node_size=20, 
    alpha=0.6
)
plt.show()

print("Number of nodes = ", G.number_of_nodes()) # 6347 nodes
print("Number of edges = ", G.number_of_edges()) # 34174 edges
print(f"Graph density: {nx.density(G)}") # 0.00169, sparse PPI network, which is expected for real biological networks. 

# --- Extract features 
# Degree - Essential genes tend to be hubs with many interactors
nx.set_node_attributes(G, dict(nx.degree(G)), 'Degree')

# Degree centrality - Normalised version of degree; scales for different graph sizes
nx.set_node_attributes(G, nx.degree_centrality(G), 'Degree Centrality')

# Closeness centrality - Captures how fast info spreads from the node
nx.set_node_attributes(G, nx.closeness_centrality(G), 'Closeness Centrality')

# Betweenness centrality - Can identify bottlenecks or relay nodes, 
nx.set_node_attributes(G, nx.betweenness_centrality(G), 'Betweenness Centrality')

# Eigenvector centrality - Detects nodes connected to other important nodes
nx.set_node_attributes(G, nx.eigenvector_centrality(G, max_iter=1000), 'Eigenvector Centrality')

# Average neighbour degree - Measures how connected the neighbors of a node are
nx.set_node_attributes(G, nx.average_neighbor_degree(G), 'Average Neighbor Degree')

# Clustering coefficient - Fraction of a node’s neighbors that are also neighbors with each other
nx.set_node_attributes(G, nx.clustering(G), 'Clustering Coefficient')

# Number of triangles - How many triangles a node is part of (3-node loops)
nx.set_node_attributes(G, nx.triangles(G), 'Triangle Count')

# K-core number - The highest k-core the node is part of. Higher values indicate more tightly-knit cores
nx.set_node_attributes(G, nx.core_number(G), 'Core Number')

# Community - Captures global position of the node in the network structure
partition = community_louvain.best_partition(G)
nx.set_node_attributes(G, partition, 'Louvain Community')

# --- Organise in dataframe 
PPIn_features = pd.DataFrame.from_dict(
    dict(G.nodes(data = True)), orient = 'index')
PPIn_features.index.name = 'FBgn_ID' # Keep gene ID labels consistent
PPIn_features.head()

# --- Filter for genes of interest only 
#  Load your gene label sets to filter features by class later
EG_ids = pd.read_csv('/Users/tanishi/FYP/FlyBase_IDs/essential_genes.txt', header=None)[0].tolist()
NEG_ids = pd.read_csv('/Users/tanishi/FYP/FlyBase_IDs/nonessential_genes.txt', header=None)[0].tolist()

# Filter EGs and NEGs separately 
PPIn_EG = PPIn_features.loc[PPIn_features.index.intersection(EG_ids)]
PPIn_NEG = PPIn_features.loc[PPIn_features.index.intersection(NEG_ids)]

# --- Save to local computer
PPIn_features.to_csv("~/FYP/feature_extraction/PPIn_all.csv")
PPIn_EG.to_csv("~/FYP/feature_extraction/PPIn_essential.csv")
PPIn_NEG.to_csv("~/FYP/feature_extraction/PPIn_nonessential.csv")

# --- Visualise to see if everything is working fine
gene_id = 'FBgn0002284'
degree_centrality_value = PPIn_features_filtered.loc[gene_id, 'Degree Centrality']
print(f"Degree centrality of {gene_id}: {degree_centrality_value}")

# Visualise degree centrality distribution across 29 genes - left skewed distribution with most genes having DC < 0.08. 3 genes have 0.14 < DC < 0.16
plt.figure(figsize=(8,5))
plt.hist(PPIn_features_filtered['Degree Centrality'], bins=20, alpha=0.7)
plt.axvline(degree_centrality_value, color='red', linestyle='dashed', linewidth=2,
            label=f'{gene_id}')
plt.xlabel('Degree Centrality')
plt.ylabel('Frequency')
plt.title('Degree Centrality Distribution of Genes')
plt.legend()
plt.show()

# Check gene's local neighbourhood in PPI graph - works!!!!!
ego = nx.ego_graph(G, gene_id)  # subgraph centered on gene_id and neighbors

plt.figure(figsize=(6,6))
pos = nx.spring_layout(ego)
nx.draw(ego, pos, with_labels=True, node_color='skyblue', edge_color='gray', node_size=500)
plt.title(f'Ego Network for {gene_id}')
plt.show()

# --- Save features as .csv
PPIn_features_filtered.to_csv("~/FYP/feature_extraction/LR/LR_PPIn_features.csv", sep = ',')