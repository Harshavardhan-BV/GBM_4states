library(trqwe)
library(ggplot2)
library(data.table)
library(dplyr)

# Read the GSE ID and sc from command line
args = commandArgs(trailingOnly=TRUE)
# test if there is only one argument: if not, return an error
if (length(args)<1) {
  stop("Atleast one argument must be supplied (GSE ID)", call.=FALSE)
}
sc = '-sc' %in% args
GSE = args[!args %in% '-sc'][1]

corr_df = function(GSE, GSM, sc){
    print(GSM)
    # Read the counts
    df = mcreadRDS(paste0("../Data_generated/", GSE, "/Counts/", GSM, "_counts.rds"), mc.cores=4) 
    # GSEA scores
    if (sc){
        df_gsea = read.csv(paste0("../Output/", GSE, "/AUCell/", GSM, "-AUCell.csv"), row.names=1)
    }else{
        df_gsea = read.csv(paste0("../Output/", GSE, "/ssGSEA/", GSM, "-ssGSEA.csv"), row.names=1)
    }
    # Select columns starting with Nef or Ver
    df_gsea = df_gsea[,grep('Nef|Ver',colnames(df_gsea))]
    # Remove rows with zero variance
    df = df[apply(df, 1, sd) != 0, ]
    # Get the correlation between AUCell and gene expression
    corrdf = data.frame(cor(t(df),df_gsea,method="spearman"))
    # Write the output
    fwrite(corrdf, paste0("../Output/",GSE,"/GRN/",GSM,"_GG-corr.tsv"), sep='\t', col.names=TRUE, row.names=TRUE)
}


# GSM file list
GSMs = list.files(paste0('../Data_generated/',GSE,'/Counts/'),'*_counts.rds') %>% gsub('_counts.rds','',.)

# Create directory for output
dir.create(paste0("../Output/", GSE, "/GRN"), showWarnings = F, recursive = T)
dir.create(paste0("../figures/", GSE, "/GRN"), showWarnings = F, recursive = T)

# Iterate over GSM samples and generate rds
for (GSM in GSMs){
    corr_df(GSE, GSM, sc)
}