import pandas as pd
import mygene
import numpy as np

# 1. Load the normalized 16-sample matrix
df = pd.read_csv('transcriptomics/tcell_exp/results/normalized_16_samples.csv', index_col=0)

# 2. Map Ensembl IDs to Gene Symbols
mg = mygene.MyGeneInfo()
ensembl_ids = df.index.tolist()

print("Querying MyGene for symbol mapping...")
# We use scopes='ensembl.gene' to match your index
results = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human')

# Create a mapping dictionary
mapping = {res['query']: res['symbol'] for res in results if 'symbol' in res}

# Add symbols to the dataframe
df['symbol'] = df.index.map(mapping)
# Remove genes that didn't map to a symbol
df_mapped = df.dropna(subset=['symbol'])

# 3. Define the T-Cell Marker Panel
# General (CD3), Cytotoxic (CD8), Helper (CD4), Regulatory (FOXP3), and Activation (GZMB, IFNG)
tcell_markers = ['CD3D', 'CD3E', 'CD8A', 'CD8B', 'CD4', 'FOXP3', 'GZMB', 'PRF1', 'IFNG', 'PDCD1']

# 4. Extract these markers
tcell_df = df_mapped[df_mapped['symbol'].isin(tcell_markers)]

# Group by symbol (in case of multiple Ensembl IDs for one gene) and take the mean
tcell_comparison = tcell_df.groupby('symbol').mean()

# 5. Save the T-cell specific profile
tcell_comparison.to_csv('transcriptomics/tcell_exp/results/tcell_marker_comparison_16.csv')

print("T-cell marker extraction complete!")
print(tcell_comparison)
