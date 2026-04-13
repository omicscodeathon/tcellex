import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

#LOAD RAW DATA (DGE requires raw integers)
# We re-combine the raw audited matrices we made earlier
luad_raw = pd.read_csv('transcriptomics/tcell_exp/results/luad_audited_matrix.csv', index_col=0)
lusc_raw = pd.read_csv('transcriptomics/tcell_exp/results/lusc_audited_matrix.csv', index_col=0)

# Merge into one 16-sample raw count matrix
counts = luad_raw.join(lusc_raw, how='inner').T # Transpose: Samples as rows
counts = counts.astype(int)

#CREATE METADATA
# This tells PyDESeq2 which sample belongs to which group
metadata = pd.DataFrame({
    'histology': ['LUAD']*8 + ['LUSC']*8
}, index=counts.index)

#RUN PyDESeq2
print("Initializing DESeq2 analysis...")
dds = DeseqDataSet(
    counts=counts, 
    metadata=metadata, 
    design="~histology" # Formula comparing the two types
)
dds.deseq2()

# EXTRACT RESULTS
# Contrast: factor, numerator, denominator
# This will show genes higher in LUSC relative to LUAD
stat_res = DeseqStats(dds, contrast=["histology", "LUSC", "LUAD"])
stat_res.summary() # Forces statistical computation
results_df = stat_res.results_df

#SAVE RESULTS
results_df.to_csv('transcriptomics/tcell_exp/results/dge_luad_vs_lusc_full.csv')

# Filter for the most significant genes (FDR < 0.05 and 2-fold change)
sig_genes = results_df[(results_df['padj'] < 0.05) & (abs(results_df['log2FoldChange']) > 1)]
sig_genes.to_csv('transcriptomics/tcell_exp/results/dge_significant_genes.csv')

print(f"Analysis complete! Found {len(sig_genes)} significant genes.")
