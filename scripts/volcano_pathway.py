import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
import mygene

#LOAD AND CLEAN DGE RESULTS
df = pd.read_csv('transcriptomics/tcell_exp/results/dge_luad_vs_lusc_full.csv', index_col=0)

# CRITICAL: Remove STAR summary rows (N_unmapped, etc.)
df = df[~df.index.str.startswith('N_')]

#MAP ENSG TO SYMBOLS
mg = mygene.MyGeneInfo()
print("Mapping Ensembl IDs to Symbols...")
map_res = mg.querymany(df.index.tolist(), scopes='ensembl.gene', fields='symbol', species='human')
mapping = {res['query']: res['symbol'] for res in map_res if 'symbol' in res}
df['symbol'] = df.index.map(mapping)

#VOLCANO PLOT
plt.figure(figsize=(10, 7))
df['neg_log10_padj'] = -np.log10(df['padj'].fillna(1))
lusc_mask = (df.pvalue < 0.05) & (df.log2FoldChange > 0.5)
luad_mask = (df.pvalue < 0.05) & (df.log2FoldChange < -0.5)

plt.scatter(df.log2FoldChange, df.neg_log10_padj, c='grey', alpha=0.3, s=10)
plt.scatter(df[lusc_mask].log2FoldChange, df[lusc_mask].neg_log10_padj, c='red', s=15, label='Up in LUSC')
plt.scatter(df[luad_mask].log2FoldChange, df[luad_mask].neg_log10_padj, c='blue', s=15, label='Up in LUAD')
plt.title("Volcano Plot: LUSC vs LUAD")
plt.savefig('transcriptomics/tcell_exp/figures/volcano_luad_vs_lusc.png')
print("Volcano plot saved.")

#PATHWAY ENRICHMENT (Using Only Mapped Symbols)
# We use cutoff=1 to force results even if significance is low
for group, condition in [("LUSC", lusc_mask), ("LUAD", luad_mask)]:
    # ONLY take valid symbols for Enrichr
    genes = df[condition]['symbol'].dropna().unique().tolist()
    
    if len(genes) > 5:
        print(f"Running Enrichment for {group} with {len(genes)} symbols...")
        enr = gp.enrichr(gene_list=genes, 
                         gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2021'], 
                         outdir=f'transcriptomics/tcell_exp/results/enr_{group}',
                         no_plot=True, 
                         cutoff=1) 
        
        if not enr.results.empty:
            # Manual barplot to ensure file creation
            res = enr.results.head(10)
            plt.figure(figsize=(12, 6))
            plt.barh(res.Term, -np.log10(res['P-value']), color='skyblue')
            plt.title(f"Top 10 Pathways in {group}")
            plt.xlabel("-log10(P-value)")
            plt.gca().invert_yaxis()
            plt.tight_layout()
            plt.savefig(f'transcriptomics/tcell_exp/figures/pathways_{group}_manual.png')
            plt.close()
            print(f"Success: {group} pathway figure saved.")
