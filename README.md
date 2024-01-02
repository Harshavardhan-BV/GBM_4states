# GBM_4states
Analysis of transcriptomic data to understand the relationship of the 4 states given in Neftel et al. and Verhaak et al. 


## Requirements
- R
    - tidyverse
    - trqwe
    - data.table
    - AUCell
    - Seurat
    - hdf5r
    - biomaRt
    - survival
    - survminer
    - ggforestplot
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
### 0. Download data
- Download the raw data from the GEO database 
- If it has a combined GSEID_RAW.tar file, extract and place the GSEID_RAW folder to the Data folder
- If the dataset has individual files, create a GSEID_RAW folder and place the files in it

### 1. Data pre-processing
- Convert to rds
Convert raw files to rds for easier loading. Done on a case-by-case basis as formats are not standardized.
```bash
cd preprocessing
Rscript counts_to_rds-GSEID.R
```

### 2. Gene-set enrichment analysis
- Get the Gene-set enrichment score for each cell/sample using AUCell/ssGSEA
```bash
cd GSEA
Rscript AUCell.R GSEID # if single cell data
#(or)
python ssGSEA.py GSEID # if bulk data
```
- Calculate the pairwise spearman correlation of the AUCell/ssGSEA scores across the cells/samples
```bash
python Corr_GSEA.py GSEID -sc #if single cell data
#(or)
python Corr_GSEA.py GSEID #if bulk data
```
#### 2.1 Heatmap of GSEA correlation (Figure 2A)
- Plot the correlation matrix as clustermap for a particular dataset
```bash
python hmap_GSEA.py GSEID
```
#### 2.2 Scatterplot of GSEA scores (Figure 2B / S1)
- Plot a scatter between each pair of subtype scores for a particular dataset 
```bash
python 2D_GSEA.py GSEID -sc #if single cell data
#(or)
python 2D_GSEA.py GSEID #if bulk data
```
#### 2.3 2D representation of GSEA scores (Figure 1A)
- Plot a Neftel style 2D representation (Ninjastar) of the AUCell scores for a particular dataset (Only for single cell data)
```bash
python Ninjastar_GSEA.py GSEID 
```
#### 2.4 Meta-analysis of GSEA correlation (Figure 4/5/S3)
- Edit the GSEs in `Input/Datsets_Bulk.csv` and `Input/Datsets_SC.csv` as necessary.
- Get the volcano plot of the correlation between NPC-MES/PN-MES 
```bash
python volc_GSEA.py
```
- Get the counts of -ve/+ve/NA correlations as a stacked barplot
```bash
python barplot_GSEA.py
```
- Get the heatmap of the correlation dataset-wise
```bash
python hmap_meta.py
```

### 3. Correlation of GBM signature expression
- Get the correlation matrix of gene expression of the signatures
```bash
cd Expression_Corr
Rscript Corr_Expr.R GSEID 
```
- Get the J-metric of correlation between each pair of subtype
```bash
python JMetric.py GSEID
```
#### 3.1 Heatmap of signature expression correlation (Figure 3A/S2A)
- Plot the correlation matrix as clustermap for a particular dataset
```bash
python hmap_Expr.py GSEID
```
#### 3.2 Barplot of J-Metric (Figure 3B/S2B)
- Plot the J-metric as a barplot for a particular dataset
```bash
python barplot_Jmetric.py GSEID
```

### 4. PCA of GBM signature expression
- Get the PCA loadings and explained variance of the GBM signature expression
```bash
cd PCA
Rscript PCA_States.R GSEID
```
#### 4.1 Plot PCA loadings (Figure 3C/S2C)
- Plot the PCA loadings and explained variance of the GBM signature expression for a particular dataset
```bash
python loadingsplot_PCA_states.py GSEID
```

### 5. TF associations 
- Get the correlation matrix of gene expression with each of the signature scores
```bash
cd GRN
Rscript Gene_GSEA_Corr.R GSEID -sc #if single cell data
#(or)
Rscript Gene_GSEA_Corr.R GSEID #if bulk data
```
#### 5.1 Plot correlation of TF with NPC/PN vs MES (Figure 6A/S4A)
- Plot the scatter plot of each gene/TF with x coordinate as correlation of TF with NPC/PN and y coordinate as correlation of TF with MES. Also saves a list of top correlated TFs
```bash
python plot_Gene_GSEA_Corr.py GSEID 
```
#### 5.2 Heatmap of top TF correlation (Figure 6B/S4B)
- Get the pairwise spearman correlation of the top correlated TFs across the cells/samples
```bash
Rscript TF_Corr.R GSEID
```
- Plot the heatmap of the top correlated TFs for a particular dataset
```bash
python plot_TF_Corr.py GSEID
```
#### 5.3 Survival analysis of top TFs (Figure 6C/S4C)
- Perform survival analysis and get the forest plot of HR and KM plots for the top correlated TFs for a particular dataset
```bash
Rscript HR_KM_TF.R GSEID
```

## Additional Miscellanous analysis
These analysis/results are not reported in the manuscript but were done for exploratory purposes. The interpretation of the results are left to the user. 

### 2. Gene-set enrichment analysis
```bash
cd GSEA
```
#### 2.5 Correlation with GSEA with gene expression
- Calculate the spearman correlation of the AUCell/ssGSEA scores with specific genes across the cells/samples
```bash
Rscript Corr_GSEA_Expr.R GSEID
```
- Plotting script not included
#### 2.6 Meta-analysis of other signatures
- Edit the GSEs in `Input/Datsets_Bulk.csv` and `Input/Datsets_SC.csv` as necessary.
- Get the volcano plot grid of the correlation for different combinations of signatures
```bash
python volc_grid_GSEA.py
```

#### 2.7 Correlation of GSEA of Extreme vs Intermediate states
- Seperate the samples into extremes: (NPC & MES or PN & MES) and intermediates (OPC & AC or NL & CL), perform correlation of AUCell/ssGSEA scores seperately for them and plots the clustermap
```bash
python Corr_split_GSEA.py GSEID
```

### 3. Correlation of GBM signature expression
```bash 
cd Expression_Corr
```
#### 3.3 Pairwise correlation
- Plot a heatmap of the correlation of the signature expression for genes in each pair of subtype for a particular dataset
```bash
python hmap_Expr_pair.py GSEID
```
#### 3.4 Meta-analysis of J-metric
- Edit the GSEs in `Input/Datsets_Bulk.csv` and `Input/Datsets_SC.csv` as necessary.
- Plot the distribution of J-metric for each pair of subtype across datasets. 
```bash
python box_Jmetric.py 
```

### 4. PCA of GBM signature expression
```bash
cd PCA
```
#### 4.2 Clustermap of PC1 and PC2 loadings
- Plot the clustermap of PC1 and PC2 loadings together for a particular dataset
```bash
python clustermap_PCA_states.py GSEID
```
#### 4.3 Pairwise PCA 
- Plot the loadings and explained variance of the PCA of the signature expression for genes in each pair of subtype for a particular dataset
```bash
python loadingsplot_PCA_pair.py GSEID
```
#### 4.4 PCA swaps
- Swap the genes in signature with housekeeping genes and get the explained variance for % of genes swapped
```bash
Rscript PCA_States_swap.R GSEID
```
- Plot the explained variance for % of genes swapped for a particular dataset
```bash
python varplot_PCA_Swaps.py GSEID
```
#### 4.5 PCA correlatin with GSEA
- Get the PCA loadings and explained variance on the expression levels of all genes
```bash
Rscript PCA_Full.R GSEID
```
- Select the top 50 +ve and -ve loadings and correlate them with AUCell/ssGSEA scores
```bash
Rscript PCA_AUCell_corr.R GSEID
```
- Plot the correlation as a clustermap
```bash
python python hmap_PCA_AUCell_corr.py GSEID
```
#### 4.6 PCA of correlation 
- Get the PCA loadings and explained variance of the correlation matrix of the signature expression
```bash
Rscript PCA_Corr.R GSEID
```
- Plot the PCA loadings and explained variance of the correlation matrix of the signature expression for a particular dataset
```bash
python loadingsplot_PCA_Corr.py GSEID
```
#### 4.7 Non-linear dimensionality reduction of GBM signature expression
- Get the 2D t-SNE and UMAP plots based on expression of GBM signature genes for a particular dataset
```bash
python tSNE_states.R GSEID
```