library(Seurat)
library(trqwe)
library(dplyr)
library(data.table)

gbm_pca = function(GSE, GSM){
    print(GSM)
    # Read the counts
    counts = mcreadRDS(paste0("./Data_generated/Imputed/", GSE, "/", GSM, "_imputed.rds"), mc.cores=4)
    if (!all(dim(counts)) || is.null(dim(counts))){return()}
    counts = counts[!apply(counts,1,sd)==0,]
    # Run PCA
    pca = prcomp(t(counts), center = T, scale. = T)
    # Get the loadings and explained variance of the first 2 PCs
    loadings = data.frame(pca$rotation[,1:2])
    exp_var = data.frame(pca$sdev[1:2]^2 / sum(pca$sdev^2))
    # Save the loadings as a tsv
    fwrite(exp_var, paste0("./Output/", GSE,'/PCA-full/',GSM,'_expvar.tsv'), sep='\t', col.names=FALSE)
    fwrite(loadings, paste0("./Output/", GSE,'/PCA-full/',GSM,'_loadings.tsv'), sep='\t', col.names=TRUE, row.names=TRUE)
}

# GSE ID
GSE = "GSE168004"

# GSM file list
GSMs = list.files(paste0("./Data_generated/", GSE, "/Imputed/"), "*_imputed.rds")  %>% gsub('_imputed.rds','',.)

# Create directory for output
dir.create(paste0("./Output/", GSE, "/PCA-full"), showWarnings = F, recursive = T)

# Iterate over GSM samples and generate rds
for (i in 1:length(GSMs)){
    gbm_pca(GSE, GSMs[i])
}
