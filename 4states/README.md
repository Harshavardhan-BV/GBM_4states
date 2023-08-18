# GBM4state
Analysis of transcriptomic data to understand the relationship of the 4 distinct states given in Neftel et al.


## Requirements
- R
    - tidyverse
    - trqwe
    - data.table
    - AUCell
    - RMagic
    - Seurat
    - hdf5r
    - biomaRt
    - 
- Python 
    - [requirements.txt](./requirements.txt)

## Structure of the repo
- [Data](./Data): Contains the raw data files (has to be downloaded seperately)
- [Data_generated](./Data_generated): Contains the pre-processed data for further analysis
- [Output](./Output): Contains the output files generated from the analysis
- [figures](./figures): Contains the figures generated from the analysis
- [preprocess](./preprocess): Contains the scripts for pre-processing the raw data in Data
- [Signature](./Signature): Contains the signature files and gene annotations

## Usage
### 1. Data pre-processing
- Convert to rds
Convert raw files to rds for easier loading. Done on a case-by-case basis as formats are not standardized.
```bash
Rscript counts_to_rds-GSEID.R
```
- Imputation (Only for single cell data)
Generate an imputed counts matrix using MAGIC.
```bash
Rscript Imputation.R GSEID
```

### 2. Gene-set enrichment analysis
- Get the Gene-set enrichment score for each cell/sample using AUCell/ssGSEA
```bash
Rscript AUCell.R GSEID # if single cell data
#(or)
python ssGSEA.py GSEID # if bulk data
```
- Calculate the correlation of the AUCell/ssGSEA scores between the samples
```bash
python Corr_GSEA.py GSEID -sc #if single cell data
#(or)
python Corr_GSEA.py GSEID #if bulk data
```
- (Optional) Plot the correlation matrix as clustermap
```bash
python hmap_GSEA.py GSEID
```
- Get the volcano plot of the correlation between combinations across samples. Edit the GSEs in `Input/Datsets_Bulk.csv` and `Input/Datsets_SC.csv`
```bash
python volcano_GSEA.py
```

### 3. Gene signature overlap
Get the overlap of the genes in each gene set
```bash
python Signature_overlap.py
```

### 4. GBM signature expression
- Get the correlation matrix of gene expression of the signatures
```bash
Rscript Corr_Expr.R GSEID -sc #if single cell data
#(or)
Rscript Corr_Expr.R GSEID #if bulk data
```
- Get the consistency of correlation between the signature expression across samples
```bash
python Consistency_Expr.py GSEID
```
- (Optional) Plot the correlation matrix as clustermap, pairwise as a heatmap, and the consistency as a barplot
```bash
python hmap_Expr.py GSEID
```
- Get the distribution of consistency for combinations. Edit the GSEs in `Input/Datsets_Bulk.csv` and `Input/Datsets_SC.csv`
```bash
python kde_Cons_Expr.py
```

### 5. PCA of GBM signature expression
Get the PCA loadings and explained variance of the GBM signature expression
```bash
Rscript PCA_States.R
python plot_PCA_states.py
```

### 6. PCA correlation with AUCell/ssGSEA scores
Get the PCA loadings and explained variance of the full expression and correlate the top/last loadings with AUCell/ssGSEA scores
```bash
Rscript PCA_Full.R
Rscript PCA_AUCell_corr.R
python plot_PCA_AUCell_corr.py
```