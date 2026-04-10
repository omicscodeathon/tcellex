import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# 1. LOAD DATA
# We use the T-cell marker comparison file created in Stage 3
tcell_exp = pd.read_csv('transcriptomics/tcell_exp/results/tcell_marker_comparison_16.csv', index_col=0)

# 2. GENERATE T-CELL INFILTRATION SCORING (BAR CHART)
# Calculate mean expression of all markers per sample
tcell_scores = tcell_exp.mean(axis=0)
score_df = pd.DataFrame({'Score': tcell_scores, 'Histology': ['LUAD']*8 + ['LUSC']*8})
score_df = score_df.sort_values('Score')

plt.figure(figsize=(12, 6))
colors = {'LUAD': 'steelblue', 'LUSC': 'darkorange'}
sns.barplot(data=score_df, x=score_df.index, y='Score', hue='Histology', palette=colors, dodge=False)
plt.xticks(rotation=90)
plt.title("T-cell Infiltration Score per Sample (LUAD vs LUSC)")
plt.ylabel("Mean log2(CPM + 1)")
plt.tight_layout()
plt.savefig('transcriptomics/tcell_exp/figures/tcell_infiltration_scores.png')
plt.close()

# 3. GENERATE COMPARATIVE HEATMAP
plt.figure(figsize=(14, 8))
# Cluster samples to see if LUAD and LUSC group separately based ONLY on T-cell markers
sns.clustermap(tcell_exp, cmap='YlGnBu', annot=True, figsize=(15, 10), 
               col_colors=['steelblue']*8 + ['darkorange']*8, # Blue for LUAD, Orange for LUSC
               standard_scale=1) # Scale genes for better visual contrast
plt.title("Heatmap: T-cell Marker Expression Profile")
plt.savefig('transcriptomics/tcell_exp/figures/tcell_heatmap.png')
plt.close()

# 4. SUBTYPE-SPECIFIC INFILTRATION BOXPLOT
plt.figure(figsize=(8, 6))
sns.boxplot(data=score_df, x='Histology', y='Score', palette=colors)
sns.stripplot(data=score_df, x='Histology', y='Score', color='black', alpha=0.5)
plt.title("Comparison of T-cell Infiltration: LUAD vs LUSC")
plt.ylabel("Immune Score (Mean Expression)")
plt.tight_layout()
plt.savefig('transcriptomics/tcell_exp/figures/tcell_subtype_comparison_boxplot.png')
plt.close()

print("Immune profiling figures generated: Scores Bar Chart, Heatmap, and Subtype Boxplot.")
