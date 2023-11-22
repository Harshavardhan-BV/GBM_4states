library(dplyr)
library(trqwe)
library(data.table)

# Read the GSE ID and sc from command line
args = commandArgs(trailingOnly=TRUE)
# test if there is only one argument: if not, return an error
if (length(args)<1) {
  stop("Atleast one argument must be supplied (GSE ID)", call.=FALSE)
}
sc = '-sc' %in% args
GSE = args[!args %in% '-sc'][1]

corr_df = function(GSE, GSM, sc, genes){
    print(GSM)
    # Read the counts
    df = mcreadRDS(paste0("../Data_generated/", GSE, "/Counts/", GSM, "_counts.rds"), mc.cores=4) 
    # GSEA scores
    if (sc){
        df_gsea = read.csv(paste0("../Output/", GSE, "/AUCell/", GSM, "-AUCell.csv"), row.names=1)
    }else{
        df_gsea = read.csv(paste0("../Output/", GSE, "/ssGSEA/", GSM, "-ssGSEA.csv"), row.names=1)
    }
    # remove NA values and those genes not in df
    genes = genes[genes %in% rownames(df)]
    # Select the rows in df based on values in gbmgenes
    df = df[genes,]

    # Get the correlation between AUCell and gene expression
    corrdf = data.frame(cor(t(df),df_gsea,method="spearman"))
    # Get the p-value of the correlation
    pvaldf = data.frame(matrix(NA, nrow=nrow(df), ncol=ncol(df_gsea)),row.names = rownames(df))
    colnames(pvaldf) = colnames(df_gsea)
    for (i in 1:nrow(df)) {
        for (j in 1:ncol(df_gsea)) {
            pvaldf[i,j] = cor.test(t(df)[,i], df_gsea[,j], method="spearman")$p.value
        }
    }
    
    # Write the output
    fwrite(corrdf, paste0("../Output/",GSE,"/GSEA-Expr/",GSM,"-corr.tsv"), sep='\t', col.names=TRUE, row.names=TRUE)
    fwrite(pvaldf, paste0("../Output/",GSE,"/GSEA-Expr/",GSM,"-pval.tsv"), sep='\t', col.names=TRUE, row.names=TRUE)
}

# GSM file list
GSMs = list.files(paste0("../Data_generated/", GSE, "/Counts/"), "*_counts.rds")  %>% gsub('_counts.rds','',.)

genes = c('CTLA4','CD274','LAG3','HAVCR2','CD47','LGALS9','CD276')
# Create directory for output
dir.create(paste0("../Output/", GSE, "/GSEA-Expr"), showWarnings = F, recursive = T)

# Iterate over GSM samples and generate rds
for (GSM in GSMs){
    corr_df(GSE, GSM, sc, genes)
}
