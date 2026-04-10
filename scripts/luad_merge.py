import pandas as pd
import os
import glob

# 1. Define paths (Update this to your LUAD data folder)
luad_dir = 'transcriptomics/tcell_exp/data/luad'
output_dir = 'transcriptomics/tcell_exp/results'

# 2. Get list of all .tsv files in the LUAD folder
luad_files = glob.glob(os.path.join(luad_dir, "*.tsv"))

if len(luad_files) < 8:
    print(f"Warning: Only {len(luad_files)} files found. Expected 8.")

# 3. Load and Merge with Audit
dfs = []
for file in luad_files:
    # TCGA STAR files typically have metadata in the first 6 rows
    df = pd.read_csv(file, sep='\t', comment='#')
    
    # Audit: Check for Null values in this specific file
    null_count = df.isnull().sum().sum()
    if null_count > 0:
        print(f"File {os.path.basename(file)} contains {null_count} null values.")
    
    # Prepare for merging: Use gene_id as index and extract 'unstranded' counts
    # Rename column to sample UUID for identification
    sample_id = os.path.basename(file).split('.')[0]
    df = df[['gene_id', 'unstranded']].set_index('gene_id')
    df.columns = [sample_id]
    dfs.append(df)

# 4. Combine all 8 LUAD samples
# We use an 'outer' join initially to catch all gene entries
luad_matrix = pd.concat(dfs, axis=1, join='outer')

# 5. Post-Merge Global Audit
print("\n--- LUAD Global Audit Summary ---")
print(f"Total Genes: {len(luad_matrix)}")
print(f"Global Null Values found: {luad_matrix.isnull().sum().sum()}")

# 6. Cleaning Step: Remove genes that are zero across all 8 samples
# This reduces noise before normalization and DGE
luad_matrix = luad_matrix[(luad_matrix.T != 0).any()]
print(f"Genes remaining after filtering zero-count rows: {len(luad_matrix)}")

# 7. Save the audited LUAD matrix
luad_matrix.to_csv(os.path.join(output_dir, 'luad_audited_matrix.csv'))
print(f"\nAudited LUAD matrix saved to {output_dir}/luad_audited_matrix.csv")
