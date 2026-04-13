import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg") # Fix for the Qt/wayland error
import matplotlib.pyplot as plt
import seaborn as sns
import mygene

#LOAD DATA
df = pd.read_csv('transcriptomics/tcell_exp/results/normalized_16_samples.csv', index_col=0)

#MAP ENSG TO SYMBOLS (Required for signature matching)
mg = mygene.MyGeneInfo()
print("Mapping Ensembl IDs to Symbols for Deconvolution...")
map_res = mg.querymany(df.index.tolist(), scopes='ensembl.gene', fields='symbol', species='human')
mapping = {res['query']: res['symbol'] for res in map_res if 'symbol' in res}
df.index = df.index.map(mapping)
df = df[df.index.notnull()] # Remove genes that didn't map

#DEFINE IMMUNE SIGNATURES
signatures = {
    'T_cells': ['CD3D', 'CD3E', 'CD2'],
    'CD8_T_cells': ['CD8A', 'CD8B'],
    'Cytotoxic': ['GZMB', 'GZMA', 'PRF1'],
    'NK_cells': ['NCAM1', 'KLRB1'],
    'B_cells': ['CD19', 'MS4A1'],
    'Macrophages': ['CD68', 'CD163']
}

#CALCULATE SCORES
decon_results = {}
for cell, markers in signatures.items():

    present_markers = [m for m in markers if m in df.index]
    if present_markers:
        print(f"Found {len(present_markers)} markers for {cell}")
        decon_results[cell] = df.loc[present_markers].mean()
    else:
        print(f"Warning: No markers found for {cell}")

if not decon_results:
    print("Error: No signature markers found in the dataset. Check mapping.")
    exit()

decon_df = pd.DataFrame(decon_results).T

#PLOT HEATMAP
plt.figure(figsize=(12, 8))
sns.heatmap(decon_df, annot=True, cmap="YlGnBu", fmt=".2f")
plt.title("Immune Cell Deconvolution: LUAD vs LUSC")
plt.ylabel("Cell Type")
plt.xlabel("Samples")
plt.tight_layout()
plt.savefig('transcriptomics/tcell_exp/figures/deconvolution_heatmap.png')
print("Success! Deconvolution heatmap saved.")
