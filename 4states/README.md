# GBM4state
Analysis of transcriptomic data to understand the relationship of the 4 distinct states given in Neftel et al.


## Requirements
- R
    - tidyverse
    - trqwe
    - data.table
    - AUCell
    - RMagic
- Python 
    - [requirements.txt](./requirements.txt)

## Usage
### 1. Data processing
Convert tsv file to rds for easier loading. Currently not an automated pipeline
```bash
Rscript counts_to_rds-GSEID.R
```
### 2. Imputation (Only for single cell data)
Generate an imputed counts matrix using MAGIC.
```bash
Rscript Imputation.R
```
### 3. Correlation of GBM signature expression
Get the correlation matrix of gene expression of the signatures
```bash
Rscript Correlation.R
python Corrplot.py
```
### 4. AUCell scoring (for single cell data)
Get the Gene enrichment score for each cell using AUCell
```bash
Rscript AUCell.R
python AUCell_plot.py
```

### 4. ssGSEA scoring (for bulk data)
Get the Gene enrichment score for each cell using ssGSEA
```bash
python ssGSEA.py
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