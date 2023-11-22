library(trqwe)
library(dplyr)
library(data.table)

# Read the GSE ID from command line
GSE = commandArgs(trailingOnly=TRUE)
# test if there is only one argument: if not, return an error
if (length(GSE)!=1) {
  stop("One argument must be supplied (GSE ID)", call.=FALSE)
} 

corr_pca = function(GSE, GSM){
    print(GSM)
    # Read the correlation matrix
    corrdf = read.delim(paste0("../Output/", GSE,'/Correlation/',GSM,'_correlation.tsv'))
    # If the correlation matrix is null, return
    if (is.null(corrdf)) {
        return()}
    rownames(corrdf) = corrdf[,1]
    corrdf = corrdf[,-1]
    # Run PCA
    pca = prcomp(corrdf, center = T, scale. = T)
    # Get the loadings and explained variance of the first 2 PCs
    loadings = data.frame(pca$rotation[,1:2])
    exp_var = data.frame(pca$sdev^2 / sum(pca$sdev^2))
    # Sort the loadings by the first PC
    loadings = loadings[order(loadings[,1], decreasing = T),]
    # Save the loadings as a tsv
    fwrite(exp_var, paste0("../Output/", GSE,'/Corr-PCA/',GSM,'_expvar.tsv'), sep='\t', col.names=FALSE)
    fwrite(loadings, paste0("../Output/", GSE,'/Corr-PCA/',GSM,'_loadings.tsv'), sep='\t', col.names=TRUE, row.names=TRUE)
    
}

# GSM file list
GSMs = list.files(paste0("../Output/", GSE, "/Correlation"), "*_correlation.tsv") %>% gsub('_correlation.tsv','',.)

# Create directory for output
dir.create(paste0("../Output/", GSE, "/Corr-PCA"), showWarnings = F, recursive = T)

# Iterate over GSM samples and generate rds
for (i in 1:length(GSMs)){
    corr_pca(GSE, GSMs[i])
}
