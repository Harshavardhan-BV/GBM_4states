library(trqwe)
library(data.table)
library(dplyr)
library(ggplot2)

# Read the GSE ID from command line
GSE = commandArgs(trailingOnly=TRUE)
# test if there is only one argument: if not, return an error
if (length(GSE)!=1) {
  stop("One argument must be supplied (GSE ID)", call.=FALSE)
} 

corrTF = function(GSE,GSM){
  # Read the counts 
  counts = mcreadRDS(paste0("../Data_generated/", GSE, "/Counts/", GSM, "_counts.rds"), mc.cores=4)
  # Read the top TFs
  TFs = list.files(paste0('../Output/',GSE,'/GRN/'),pattern = paste0(GSM,'_.*-topTFs.tsv'), full.names = TRUE) %>% lapply(read.delim, header=F) %>% unlist() %>% unique()
  # Subset the TFs
  TFs = TFs[TFs %in% rownames(counts)]
  counts = counts[TFs,]
  # Remove cases where standard deviation is zero
  counts <- counts[ apply(counts, 1, sd) != 0, ]
  # Get the correlation each TF
  corrdf = cor(t(counts),method="spearman") %>% data.frame()
  # Save the correlation matrix
  fwrite(corrdf, paste0("../Output/", GSE,'/GRN/',GSM,'_TF-corr.tsv'), sep='\t', col.names=TRUE, row.names=TRUE)
}

# GSM file list
GSMs = list.files(paste0('../Data_generated/',GSE,'/Counts/'),'*_counts.rds') %>% gsub('_counts.rds','',.)

# Run the function for each GSM
for (GSM in GSMs){
  corrTF(GSE,GSM)
}
