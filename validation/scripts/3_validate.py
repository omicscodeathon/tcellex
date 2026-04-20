import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib
matplotlib.use('Agg') 

#SETUP MODULAR PATHS
val_dir = os.path.expanduser("~/transcriptomics/tcell_exp/validation")
data_dir = os.path.join(val_dir, "data")
fig_dir = os.path.join(val_dir, "figures")
res_dir = os.path.join(val_dir, "results")

#Ensure results/figures folders exist
os.makedirs(fig_dir, exist_ok=True)
os.makedirs(res_dir, exist_ok=True)

#Raw files locations
matrix_file = os.path.join(data_dir, "GSE127465_human_counts_normalized_54773x41861.mtx.gz")
genes_file = os.path.join(data_dir, "GSE127465_gene_names_human_41861.tsv.gz")
meta_file = os.path.join(data_dir, "GSE127465_human_cell_metadata_54773x25.tsv.gz")

#Set scanpy figure directory explicitly
sc.settings.figdir = fig_dir

print("--- STAGE 1: LOADING & INTEGRATION ---")
#Load Sparse Matrix
adata = sc.read_mtx(matrix_file)

#Shape Alignment Check
if adata.shape[0] != 54773:
    print(f"Current shape {adata.shape}. Transposing to (Cells, Genes)...")
    adata = adata.transpose()

#Attach Gene Names
genes = pd.read_csv(genes_file, header=None, sep='\t').squeeze().astype(str).tolist()
adata.var_names = genes
adata.var_names_make_unique()
print("Gene symbols assigned successfully.")

#Positional Metadata Integration

metadata = pd.read_csv(meta_file, sep='\t')
if len(metadata) == adata.n_obs:
    #Synchronize index names first to avoid AnnData conflicts
    adata.obs_names = metadata.index.astype(str)
    adata.obs = metadata
    print("Metadata integrated successfully via positional alignment.")
else:
    print(f"Error: Count mismatch. Meta: {len(metadata)}, Matrix: {adata.n_obs}")

print("--- STAGE 2: MAPPING HISTOLOGY (LUAD VS LUSC) ---")
#Mapping Patient IDs to histology labels (Based on Zilionis et al. 2019) p1-3 = LUAD, p4-7 = LUSC
histology_map = {
    'p1': 'LUAD', 'p2': 'LUAD', 'p3': 'LUAD',
    'p4': 'LUSC', 'p5': 'LUSC', 'p6': 'LUSC', 'p7': 'LUSC'
}
adata.obs['histology'] = adata.obs['Patient'].map(histology_map)

print("--- STAGE 3: CLEANING & QUALITY CONTROL ---")
sc.pp.filter_cells(adata, min_genes=200)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
#Retain high-quality cells (< 5% mitochondrial transcripts)
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

print("--- STAGE 4: NORMALIZATION & LOG-TRANSFORM ---")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

print("--- STAGE 5: DIMENSIONALITY REDUCTION ---")
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

print("--- STAGE 6: FINAL VISUALIZATION ---")
#Generates the 'Gold Standard' proof for your presentation
sc.pl.umap(adata, color=['histology', 'CD3D'], 
           title=['Tumor Histology (LUAD vs LUSC)', 'T-Cell Signal (CD3D)'],
           palette={'LUAD': '#1f77b4', 'LUSC': '#ff7f0e'}, 
           save='_Final_Validation_Atlas.png', show=False)

#Save the finalized model for results archival
h5ad_path = os.path.join(data_dir, "GSE127465_final_validated.h5ad")
adata.write(h5ad_path)

print(f"SUCCESS! Final Validation Shape: {adata.shape}")
print(f"Final Image: {fig_dir}/umap_Final_Validation_Atlas.png")
