import pandas as pd
import os
import matplotlib.pyplot as plt

#Setup paths
input_dir = "/home/hp/transcriptomics/tcell_exp/PPI/files/"
output_dir = "/home/hp/transcriptomics/tcell_exp/results/"

#10-Gene Signature
canonical_markers = ['CD3D', 'CD3E', 'CD3G', 'CD8A', 'CD8B', 'CD4', 'GZMB', 'GZMA', 'PRF1', 'IFNG']

files = {
    'Betweenness': 'between.csv', 
    'Closeness': 'close.csv', 
    'Degree': 'degree.csv', 
    'MNC': 'MNC.csv', 
    'MCC': 'MCC.csv'
}

def run_ppi_audit(input_path, output_path):
    
    master_df = pd.DataFrame({'gene': canonical_markers})
    master_df['gene'] = master_df['gene'].astype(str)
    
    print("--- Executing Type-Safe Integrated Audit ---")
    
    for metric, filename in files.items():
        file_path = os.path.join(input_path, filename)
        if not os.path.exists(file_path): 
            print(f"Skipping: {filename} (not found)")
            continue
        
        df = pd.read_csv(file_path)
        
        #Pull node_name and the score
        val_col = metric if metric in df.columns else df.columns[-1]
        subset = df[['node_name', val_col]].copy()
        subset.columns = ['gene', metric]
        
        
        subset['gene'] = subset['gene'].astype(str)
        
        #Merge only our 10 markers
        master_df = pd.merge(master_df, subset, on='gene', how='left').fillna(0)

    #Ranking logic
    metrics_found = [m for m in files.keys() if m in master_df.columns]
    for m in metrics_found:
        # Rank the best 
        master_df[f'{m}_Rank'] = master_df[m].rank(ascending=False, method='min')
    
    rank_cols = [c for c in master_df.columns if '_Rank' in c]
    #Invert the rank so higher bars
    master_df['Integrated_Hub_Strength'] = 11 - master_df[rank_cols].mean(axis=1)
    
    # Sort for plotting
    master_df = master_df.sort_values('Integrated_Hub_Strength')

    #Generating the Summary Plot
    plt.figure(figsize=(12, 7))
    bars = plt.barh(master_df['gene'], master_df['Integrated_Hub_Strength'], 
                    color='skyblue', edgecolor='navy', alpha=0.8)
    
    #Add numerical labels to the bars
    for bar in bars:
        width = bar.get_width()
        plt.text(width + 0.1, bar.get_y() + bar.get_height()/2, 
                 f'{width:.1f}', va='center', fontweight='bold')

    plt.title('Top 10 Executive Hub Genes: Integrated Centrality Audit', fontsize=14, fontweight='bold')
    plt.xlabel('Hub Strength (Composite Rank across 5 Metrics)', fontsize=12)
    plt.ylabel('T-Cell Canonical Markers', fontsize=12)
    plt.xlim(0, 12) 
    plt.grid(axis='x', linestyle='--', alpha=0.4)
    plt.tight_layout()
    
    plot_path = os.path.join(output_path, "ppi_canonical_hub_plot.png")
    plt.savefig(plot_path, dpi=300)
    print(f"\nSUCCESS: Integrated plot saved to: {plot_path}")

if __name__ == '__main__':
    run_ppi_audit(input_dir, output_dir)
