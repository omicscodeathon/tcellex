import scanpy as sc
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

#Setup Paths
val_dir = os.path.expanduser("~/transcriptomics/tcell_exp/validation")
res_dir = os.path.join(val_dir, "results")
fig_dir = os.path.join(val_dir, "figures")

#Load the finalized model
adata = sc.read_h5ad(os.path.join(val_dir, "data/GSE127465_final_validated.h5ad"))

print("--- STARTING QUANTITATIVE ABUNDANCE AUDIT ---")

#Calculate Cellular Abundance

abundance = adata.obs.groupby(['histology', 'Major cell type'], observed=True).size().unstack(fill_value=0)
abundance_pct = abundance.div(abundance.sum(axis=1), axis=0) * 100

#Identify T-cell Column and Extract Proof
tcell_col = [col for col in abundance_pct.columns if 'T cell' in col or 'Lymph' in col][0]
proof_data = abundance_pct[tcell_col].reset_index()
proof_data.columns = ['Histology', 'T_Cell_Percentage']

print(f"Detected T-cell column: {tcell_col}")
print(proof_data)

#Generate Professional Bar Plot
plt.figure(figsize=(8, 6))
sns.set_style("whitegrid")

#Define professional colors: LUAD (Blue), LUSC (Orange)
colors = {'LUAD': '#1f77b4', 'LUSC': '#ff7f0e'}

ax = sns.barplot(x='Histology', y='T_Cell_Percentage', data=proof_data, 
                 palette=colors, edgecolor='black', linewidth=1.5)

#Add data labels on top of bars
for p in ax.patches:
    ax.annotate(f'{p.get_height():.1f}%', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha = 'center', va = 'center', 
                xytext = (0, 9), 
                textcoords = 'offset points',
                fontsize=12, fontweight='bold')

plt.title("Quantitative Validation: T-Cell Infiltration Gap", fontsize=15, pad=20)
plt.ylabel("Percentage of Total Cellularity (%)", fontsize=12)
plt.xlabel("Histological Subtype", fontsize=12)
plt.ylim(0, max(proof_data['T_Cell_Percentage']) + 5) # Add space for labels

#Save Outputs
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "Tcell_Abundance_Validation_Barplot.png"), dpi=300)
proof_data.to_csv(os.path.join(res_dir, "final_abundance_proof.csv"), index=False)

print(f"\nSUCCESS!")
print(f"Bar Plot saved to: {fig_dir}/Tcell_Abundance_Validation_Barplot.png")
print(f"Data Table saved to: {res_dir}/final_abundance_proof.csv")
