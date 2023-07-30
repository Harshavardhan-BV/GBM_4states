library(Rmagic)
library(trqwe)
library(dplyr)

# Read the GSE ID from command line
GSE = commandArgs(trailingOnly=TRUE)
# test if there is only one argument: if not, return an error
if (length(GSE)!=1) {
  stop("One argument must be supplied (GSE ID)", call.=FALSE)
} 

impute <- function(GSE,GSM) {
    count_matrix = mcreadRDS(paste0("../Data_generated/",GSE,"/Counts/",GSM,"_counts.rds"), mc.cores=4)
    count_matrix = t(count_matrix)
    gc()
    count_matrix = magic(count_matrix, n.jobs = -2, solver='exact')
    gc()
    mcsaveRDS(t(count_matrix$result), paste0("../Data_generated/",GSE,"/Imputed/",GSM,"_imputed.rds"), mc.cores = 4)
}

# GSM file list
GSMs = list.files(paste0('../Data_generated/',GSE,'/Counts/'),'*_counts.rds') %>% gsub('_counts.rds','',.)
print(GSMs)
# Create directory for output
dir.create(paste0("../Data_generated/", GSE, "/Imputed"), showWarnings = F, recursive = T)

# Iterate over GSM samples and generate rds
for (GSM in GSMs){
    impute(GSE,GSM)
}

