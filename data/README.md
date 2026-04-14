DATA DIRECTORY

OVERVIEW

This directory contains RNA-seq gene expression datasets used for the analysis of T-cell expression in lung cancer
The study focuses on a comparative analysis between Lung Adenocarcinoma(LUAD) and Lung Squamous Cell carcinoma(LUSC)
All samples included are only tumor samples thus enabling subtype specific immune analysis

DATA SOURCE

Data were obtained from The Cancer genome Atlas (TCGA). This data was accessed via the Genomic Data Commons (GDC) data portal

DATA TYPE:

-RNA-seq gene expression (STAR gene counts)
-File format : "*.tsv"

DIRECTORY STRUCTURE

data/
   luad/-contains luad tsv files
   lusc/-contains lusc tsv files

FILE FORMAT DESCRIPTION

Each .tsv file contains gene-level expression data with the following
 - gene_id- Ensembl gene identifier
 - gene_name- Gene symbol
 - gene_type- gene classification
 - unstranded- raw read counts
 - stranded_first/ stranded_second -strand specific counts
 - TPM/FPKM values- normalized measures

DATA PROCESSING STEPS

Below are the preprocessing steps used:
 - Filtering: non-gene rows were removed(N_umapped, Nmultimapping, N_noFeature, N_ambiguous)
 - Column specification : the unstranded count column was retained for downstream analysis
 - Gene handling : the gene IDs were used as primary identifiers and duplicate gene entries were removed to ensure uniqueness
  
DATA ACCESS

The files were downloaded using the GDC API uuid-based links:
  https://api.gdc.cancer.gov/data/<filename>
  each filename corresponds to a unique TCGA sample
  
⚠️ NOTE:

- All samples are tumor samples
- this study supports subtype comparison and not tumor vs normal comparison
- Data represent bulk RNA-seq not single-cell sequencing
- sample size is limited to 8 per subtype

DATA ROLE IN PIPELINE ANALYSIS

- sample merging into expression matrices
- Normalization (CPM and log2 transformation)
- Differential expression analysis
- T-cell marker extraction
- Immune scoring and clustering

DATA USAGE POLICY

All data are publicly available through TCGA and are used strictly for academic and research purposes

SUMMARY:

This dataset enables:
- Comparative transcriptomics analysis between LUAD and LUSC
- Investigation of T-cell related immune activity
- characterization of tumor immune microenvironment 





