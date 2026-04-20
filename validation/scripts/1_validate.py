import scanpy as sc
import pandas as pd
import os

#Setup Paths
val_dir = os.path.expanduser("~/transcriptomics/tcell_exp/validation")
data_dir = os.path.join(val_dir, "data")

matrix_path = os.path.join(data_dir, "GSE127465_human_counts_normalized_54773x41861.mtx.gz")
genes_path = os.path.join(data_dir, "GSE127465_gene_names_human_41861.tsv.gz")
metadata_path = os.path.join(data_dir, "GSE127465_human_cell_metadata_54773x25.tsv.gz")

print("Stage 1: Building the Validation Object...")

#Load the Sparse Matrix
adata = sc.read_mtx(matrix_path)

#Shape Alignment Check (Transpose if needed)
if adata.shape[0] != 54773:
    print(f"Current shape {adata.shape}. Transposing to (Cells, Genes)...")
    adata = adata.transpose()

#Attach Gene Names
#Squeeze into a list and make them unique to prevent plotting errors
genes = pd.read_csv(genes_path, header=None, sep='\t').squeeze().astype(str).tolist()
adata.var_names = genes
adata.var_names_make_unique()
print("Gene symbols assigned and made unique.")

#Attach Cell Metadata (The Fix)
metadata = pd.read_csv(metadata_path, sep='\t')

#ALIGNMENT FIX: 
#Instead of matching on names, we assign the metadata directly because 
#we know the order of cells (54773) matches between both files.
if len(metadata) == adata.n_obs:
    #Set the AnnData index to match the metadata index
    adata.obs_names = metadata.index.astype(str)
    adata.obs = metadata
    print("Metadata successfully integrated by position.")
else:
    print(f"Error: Row mismatch. Metadata has {len(metadata)} rows, Matrix has {adata.n_obs} cells.")

#Save as h5ad
output_path = os.path.join(data_dir, "GSE127465_integrated.h5ad")
adata.write(output_path)

print(f"SUCCESS! Validation Set Ready.")
print(f"Final Shape: {adata.shape} (Cells, Genes)")
