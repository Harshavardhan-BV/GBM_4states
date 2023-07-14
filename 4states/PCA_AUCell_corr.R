library(Seurat)
library(trqwe)
library(dplyr)
library(data.table)

load_corr_specific = function(top50_expr, last50_expr, auc, suff, GSE, GSM){
    auc = auc[,grep(suff, colnames(auc))]
    # Get the correlation of t(auc) with top50_expr and last50_expr
    top50_corr =  data.frame(cor(t(top50_expr), auc, method = 'spearman'))
    last50_corr =  data.frame(cor(t(last50_expr), auc, method = 'spearman'))
    # Write the correlation to file
    fwrite(top50_corr, paste0("./Output/", GSE,'/PCA-full/',GSM,'_top50-corr_',suff,'.tsv'), sep='\t', col.names=TRUE, row.names=TRUE)
    fwrite(last50_corr, paste0("./Output/", GSE,'/PCA-full/',GSM,'_last50-corr_',suff,'.tsv'), sep='\t', col.names=TRUE, row.names=TRUE)

}

top_load_corr = function(GSE, GSM, sc){
    print(GSM)
    # Read the counts
    if (sc){
        counts = mcreadRDS(paste0("./Data_generated/", GSE, "/Imputed/", GSM, "_imputed.rds"), mc.cores=4)
    }else{
        counts = mcreadRDS(paste0("./Data_generated/", GSE, "/Counts/", GSM, "_counts.rds"), mc.cores=4)
    }
    if (!all(dim(counts)) || is.null(dim(counts))){return()}
    # Read the loadings
    loads = read.delim(paste0("./Output/", GSE,'/PCA-full/',GSM,'_loadings.tsv'))
    # sort loads according to PC1
    loads = loads[order(loads[,2], decreasing = T),]
    # Get the top 50 genes and last 50 genes
    top50 = loads[1:50,1]
    last50 = loads[(nrow(loads)-49):nrow(loads),1]
    # Get the expression of these genes
    top50_expr = counts[top50,]
    last50_expr = counts[last50,]
    # Read the AUCell scores
    if (sc){
        auc = read.csv(paste0("./Output/", GSE,'/AUCell/',GSM,'-AUCell.csv'))
    }else{
        auc = read.csv(paste0("./Output/", GSE,'/ssGSEA/',GSM,'-ssgsea.csv'))
    }
    rownames(auc) = auc[,1]
    load_corr_specific(top50_expr, last50_expr, auc, 'Nef', GSE, GSM)
    load_corr_specific(top50_expr, last50_expr, auc, 'Ver', GSE, GSM)
    
}

# GSE ID
# GSE = "GSE168004"
# GSE = "GSE131928"
# GSE = "GSE182109"
# sc = TRUE

GSE = "CCLE"
# GSE = "TCGA"
sc = FALSE

# Get all the files in folder
GSMs = list.files(paste0("./Output/", GSE, "/PCA-full/"), "*_loadings.tsv")  %>% gsub('_loadings.tsv','',.)

# Iterate over GSM samples and generate rds
for (i in 1:length(GSMs)){
    top_load_corr(GSE, GSMs[i], sc)
}