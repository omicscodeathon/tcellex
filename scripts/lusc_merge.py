import pandas as pd
import os
import glob

#Define paths (Update this to your LUSC data folder)
lusc_dir = 'transcriptomics/tcell_exp/data/lusc'
output_dir = 'transcriptomics/tcell_exp/results'

#Get list of all .tsv files in the LUAD folder
lusc_files = glob.glob(os.path.join(lusc_dir, "*.tsv"))

if len(lusc_files) < 8:
    print(f"Warning: Only {len(lusc_files)} files found. Expected 8.")

#Load and Merge with Audit
dfs = []
for file in lusc_files:
    #TCGA STAR files typically have metadata in the first 6 rows
    df = pd.read_csv(file, sep='\t', comment='#')
    
    #Audit: Check for Null values in this specific file
    null_count = df.isnull().sum().sum()
    if null_count > 0:
        print(f"File {os.path.basename(file)} contains {null_count} null values.")
    
    #Prepare for merging: Use gene_id as index and extract 'unstranded' counts
    #Rename column to sample UUID for identification
    sample_id = os.path.basename(file).split('.')[0]
    df = df[['gene_id', 'unstranded']].set_index('gene_id')
    df.columns = [sample_id]
    dfs.append(df)

#Combine all 8 LUSC samples, We use an 'outer' join initially to catch all gene entries
lusc_matrix = pd.concat(dfs, axis=1, join='outer')

#Post-Merge Global Audit
print("\n--- LUSC Global Audit Summary ---")
print(f"Total Genes: {len(lusc_matrix)}")
print(f"Global Null Values found: {lusc_matrix.isnull().sum().sum()}")

#Cleaning Step: Remove genes that are zero across all 8 samples,This reduces noise before normalization and DGE
luad_matrix = lusc_matrix[(lusc_matrix.T != 0).any()]
print(f"Genes remaining after filtering zero-count rows: {len(luad_matrix)}")

#Save the audited LUSC matrix
lusc_matrix.to_csv(os.path.join(output_dir, 'lusc_audited_matrix.csv'))
print(f"\nAudited LUSC matrix saved to {output_dir}/lusc_audited_matrix.csv")
