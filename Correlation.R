library(dplyr)
library(trqwe)
library(data.table)
library(readxl)

corr_df = function(GSE, GSM, sc){
    print(GSM)
    # Read the counts
    if (sc){
        df = mcreadRDS(paste0("./Data_generated/", GSE, "/Imputed/", GSM, "_imputed.rds"), mc.cores=4)
    }else{
        df = mcreadRDS(paste0("./Data_generated/", GSE, "/Counts/", GSM, "_counts.rds"), mc.cores=4)
    }
    # Get the signature genes
    gbmgenes = read.csv("./Signatures/GBM_signatures.csv")
    # add values in all columns of gbmgenes to a list
    genes = unlist(gbmgenes, use.names = FALSE)
    # remove NA values and those genes not in df
    genes = genes[! is.na(genes) & (genes %in% rownames(df))]
    # Select the rows in df based on values in gbmgenes
    df = df[genes,]
    # Remove rows with 0 sd
    df = df[!apply(df,1,sd)==0,]
    # Get the correlation
    corrdf = data.frame(cor(t(df),method="spearman"))
    # Write the output
    fwrite(corrdf, paste0("./Output/",GSE,"/Correlation/",GSM,"_correlation.tsv"), sep='\t', col.names=TRUE, row.names=TRUE)
}

# GSE ID
GSE = "GSE168004"
# GSE = "GSE131928"
# GSE = "GSE182109"
sc = TRUE

# GSE = "CCLE"
# GSE = "TCGA"
# sc = FALSE

# GSM file list
if (sc){
    GSMs = list.files(paste0("./Data_generated/", GSE, "/Imputed/"), "*_imputed.rds")  %>% gsub('_imputed.rds','',.)
}else{
    GSMs = list.files(paste0("./Data_generated/", GSE, "/Counts/"), "*_counts.rds")  %>% gsub('_counts.rds','',.)
}

# Create directory for output
dir.create(paste0("./Output/", GSE, "/Correlation"), showWarnings = F, recursive = T)

# Iterate over GSM samples and generate rds
for (i in 1:length(GSMs)){
    corr_df(GSE, GSMs[i], sc)
}
