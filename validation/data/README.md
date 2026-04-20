     VALIDATION DATA(scRNA-seq)

     OVERVIEW

This directory contains data used for single-cell RNA sequencing(scRNA-seq) validation of T-cell activity in lung cancer

The data supports validation of the initial bulk RNA-seq findings by providing cell-level gene expression information

     DATASETS SOURCE

The data was obtained from the Gene Expression Omnibus(GEO):

-Acession: GSE127465

-Study: Single-cell RNA-seq of human lung tumor microenvironment

     DESCRIPTION

1. Gene names:

File: GSE127465_gene_names_human_41861.tsv.gz

Description: It contains gene identifiers corresponding to the expression matrix

2. Cell metadata:
 
File: GSE127465_human_cell_metadata_54773x25.tsv.gz

Description: Contains annotations for each cell including metadata such as cell type and sample information.

3. Expression matrix:

File: GSE127465_human_counts_normalized_54773x41861.mtx.gz

Description: Sparse matrix of normalized gene expression values

Dimensions: 54,773 cellx41,861 genes

     DATA FORMAT

-matrix format: .mtx(sparse matrix format)

-Metadata: .tsv(tab-separated values)

-cempressed format: .gz

     USAGE

These files are used as input for validation scripts locateds in the validation/scripts/ directory.

They are required to:

- load single-cell expression data

- Map gene identifiesr

- Compute T-cell signature scores

- Generate validation visualizations

     NOTES

-The dataset is large and may require subsetting for analysis on systems with limited memory

- No modifications were made to the original downloaded files

- Users are encouraged to downnload the data directly from GEO.

     ACCESS

Due to size limitations, the data might not be include in this repository. Therefore they can be downloaded from:

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127465
