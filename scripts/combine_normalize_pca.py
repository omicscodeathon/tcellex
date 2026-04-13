import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")  # Set backend for headless environments
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import os

#Load the audited matrices
luad = pd.read_csv('transcriptomics/tcell_exp/results/luad_audited_matrix.csv', index_col=0)
lusc = pd.read_csv('transcriptomics/tcell_exp/results/lusc_audited_matrix.csv', index_col=0)

#Grand Combination (Inner Join)
combined_df = luad.join(lusc, how='inner')

# Strip Ensembl version suffixes correctly
combined_df.index = combined_df.index.str.split('.').str.get(0)

#Global Normalization: log2(CPM + 1)
cpm = (combined_df / combined_df.sum()) * 1e6
log_cpm = np.log2(cpm + 1)

# Ensure figures directory exists
os.makedirs('transcriptomics/tcell_exp/figures', exist_ok=True)

#VISUALIZATION 1: BOXPLOT 
plt.figure(figsize=(14, 6))
melted_df = log_cpm.melt(var_name='Sample', value_name='Expression')
melted_df['Group'] = melted_df['Sample'].apply(lambda x: 'LUAD' if x in luad.columns else 'LUSC')

sns.boxplot(data=melted_df, x='Sample', y='Expression', hue='Group', dodge=False)
plt.xticks(rotation=90)
plt.title("Post-Normalization Expression Distribution (16 Samples)")
plt.ylabel("log2(CPM + 1)")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig('transcriptomics/tcell_exp/figures/normalization_boxplot.png')
plt.close() 

# VISUALIZATION 2: PCA PLOT 
pca = PCA(n_components=2)
components = pca.fit_transform(log_cpm.T)

pca_df = pd.DataFrame(data=components, columns=['PC1', 'PC2'], index=log_cpm.columns)
pca_df['Histology'] = ['LUAD' if x in luad.columns else 'LUSC' for x in pca_df.index]

plt.figure(figsize=(10, 8))
sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Histology', s=100, style='Histology')

for i, txt in enumerate(pca_df.index):
    plt.annotate(txt, (pca_df.PC1.iloc[i], pca_df.PC2.iloc[i]), fontsize=8, alpha=0.7)

plt.title(f"PCA: LUAD vs LUSC Clustering (PC1 Var: {pca.explained_variance_ratio_[0]:.2%})")
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.2%})")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.2%})")
plt.tight_layout()
plt.savefig('transcriptomics/tcell_exp/figures/pca_16_samples.png')
plt.close()

#Save the final matrix
log_cpm.to_csv('transcriptomics/tcell_exp/results/normalized_16_samples.csv')
print("Grand Merge Complete! Boxplot and PCA saved to transcriptomics/tcell_exp/figures/")
