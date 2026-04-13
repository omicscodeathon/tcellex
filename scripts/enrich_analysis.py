import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gseapy as gp
import mygene
import os

#LOAD DATA
results_csv = 'transcriptomics/tcell_exp/results/dge_luad_vs_lusc_full.csv'
df = pd.read_csv(results_csv, index_col=0)

#CRITICAL: Clean the IDs (ENSG00000002586.15 -> ENSG00000002586)
# This is why your previous mapping failed
df.index = df.index.astype(str).str.split('.').str.get(0)
df = df[~df.index.str.startswith('N_')] # Remove STAR noise

#MAP TO SYMBOLS
print("Mapping Cleaned Ensembl IDs to Symbols...")
mg = mygene.MyGeneInfo()
map_res = mg.querymany(df.index.tolist(), scopes='ensembl.gene', fields='symbol', species='human')
mapping = {res['query']: res['symbol'] for res in map_res if 'symbol' in res}
df['symbol'] = df.index.map(mapping)

# RUN ENRICHMENT
# Using relaxed thresholds (p < 0.05 and LFC > 0.5) to ensure we find pathways
for group, condition in [("LUSC_High", df.log2FoldChange > 0.5), 
                         ("LUAD_High", df.log2FoldChange < -0.5)]:
    
    # Filter for significant genes with a valid symbol
    gene_list = df[condition & (df.pvalue < 0.05)]['symbol'].dropna().unique().tolist()
    
    if len(gene_list) > 5:
        print(f"Running Enrichment for {group} with {len(gene_list)} genes...")
        enr = gp.enrichr(gene_list=gene_list,
                         gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2021'],
                         organism='human',
                         no_plot=True,
                         cutoff=1)
        
        if not enr.results.empty:
            # Manual Barplot to ensure the PNG is saved
            res = enr.results.sort_values('P-value').head(10)
            plt.figure(figsize=(10, 8))
            plt.barh(res.Term, -np.log10(res['P-value']), color='skyblue', edgecolor='black')
            plt.title(f"Top Pathways: {group}")
            plt.xlabel("-log10(P-value)")
            plt.gca().invert_yaxis()
            plt.tight_layout()
            
            save_path = f'transcriptomics/tcell_exp/figures/pathways_{group}_final.png'
            plt.savefig(save_path)
            plt.close()
            print(f"SUCCESS: Saved {save_path}")
    else:
        print(f"Skipping {group}: Only found {len(gene_list)} significant symbols.")
