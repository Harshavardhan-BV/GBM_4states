library(Seurat)
library(trqwe)
library(dplyr)
library(data.table)

get_sig = function(sigs,suff){
    # If suff doesnt contain - then return the genes in the signature set
    if (!grepl('-', suff)){
        sigs = sigs[,grep(suff, colnames(sigs))]
    } else {
    # Split by - and return the genes in each signature
        sigs = lapply(strsplit(suff,'-')[[1]], function(x) sigs[,grep(x, colnames(sigs))])
    }
    return(sigs)
}

gbm_pca = function(GSE, GSM, suff){
    print(paste(GSM, suff))
    # Read the counts
    counts = mcreadRDS(paste0("./Data_generated/", GSE, "/Counts/", GSM, "_counts.rds"), mc.cores=4)
    # Read the signatures
    sigs = read.csv("./Signatures/GBM_signatures.csv")
    #split suff by - and get the genes in each signature
    sigs = get_sig(sigs, suff)
    genes = unique(unlist(sigs, use.names = F))
    genes = genes[genes %in% rownames(counts)]
    # Subset only the signature genes
    counts = counts[genes,]
    if (!all(dim(counts)) || is.null(dim(counts))){return()}
    counts = counts[!apply(counts,1,sd)==0,]
    # Run PCA
    pca = prcomp(t(counts), center = T, scale. = T)
    # Get the loadings and explained variance of the first 2 PCs
    loadings = data.frame(pca$rotation[,1:2])
    exp_var = data.frame(pca$sdev^2 / sum(pca$sdev^2))
    # Save the loadings as a tsv
    fwrite(exp_var, paste0("./Output/", GSE,'/PCA/',GSM,'_expvar_',suff,'.tsv'), sep='\t', col.names=FALSE)
    fwrite(loadings, paste0("./Output/", GSE,'/PCA/',GSM,'_loadings_',suff,'.tsv'), sep='\t', col.names=TRUE, row.names=TRUE)
}

gbm_pair_pca = function(GSE, GSM, suff){
    # Read the signatures
    sigs = read.csv("./Signatures/GBM_signatures.csv")
    sigs = sigs[,grep(suff, colnames(sigs))]
    # Make combinations of the signatures
    sigs = combn(colnames(sigs),2)
    # Iterate over the combinations
    for (i in 1:ncol(sigs)){
        suff = sigs[,i]
        suff = paste(suff, collapse = '-')
        gbm_pca(GSE, GSM, suff)
    }
}

# GSE ID
# GSE = "GSE168004"
GSE = "GSE131928"
# GSE = "GSE182109"

# GSE = "CCLE"
# GSE = "TCGA"

# GSM file list
GSMs = list.files(paste0("./Data_generated/", GSE, "/Counts/"), "*_counts.rds")  %>% gsub('_counts.rds','',.)

# Create directory for output
dir.create(paste0("./Output/", GSE, "/PCA"), showWarnings = F, recursive = T)

# Iterate over GSM samples and generate rds
for (i in 1:length(GSMs)){
    gbm_pca(GSE, GSMs[i], 'Nef')
    gbm_pair_pca(GSE, GSMs[i], 'Nef')
    gbm_pca(GSE, GSMs[i], 'Ver')
    gbm_pair_pca(GSE, GSMs[i], 'Ver')
}
