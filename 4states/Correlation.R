library(dplyr)
library(trqwe)
library(data.table)
library(readxl)

corr_df = function(GSE, GSM){
    df = mcreadRDS(paste0("./Data_generated/",GSE,"/",GSM,"_imputed.rds"), mc.cores=4)
    gbmgenes = read_excel("./Signatures/1-s2.0-S0092867419306877-mmc2.xlsx",skip = 4)
    gbmgenes = gbmgenes[,0:6]
    # add values in all columns of gbmgenes to a list
    genes = unlist(gbmgenes, use.names = FALSE)
    # remove NA values and those genes not in df
    genes = genes[! is.na(genes) & (genes %in% rownames(df))]
    # Select the rows in df based on values in gbmgenes
    df = df[genes,]
    # Remove rows with 0 sd
    df = df[!apply(df,1,sd)==0,]
    corrdf = data.frame(cor(t(df),method="spearman"))
    fwrite(corrdf, paste0("./Output/",GSE,"/",GSM,"_correlation.tsv"), sep='\t', col.names=TRUE, row.names=TRUE)
}

# GSE ID
GSE = "GSE168004"

# GSM file list
GSMs = c("OSM_celllines", "mgg23", "mgg75")

# Iterate over GSM samples and generate rds
for (i in 1:length(GSMs)){
    corr_df(GSE, GSMs[i])
}
