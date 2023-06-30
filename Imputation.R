library(Rmagic)
library(trqwe)

impute <- function(GSE,GSM) {
    count_matrix = mcreadRDS(paste0("./Data_generated/",GSE,"/Counts/",GSM,"_counts.rds"), mc.cores=4)
    count_matrix = t(count_matrix)
    gc()
    count_matrix = magic(count_matrix, n.jobs = -2, solver='exact')
    gc()
    mcsaveRDS(t(count_matrix$result), paste0("./Data_generated/",GSE,"/Imputed/",GSM,"_imputed.rds"), mc.cores = 4)
}

# GSE ID
GSE = "GSE168004"

# GSM file list
GSMs = list.files(paste0('./Data_generated/',GSE,'/Counts/'),'*_counts.rds') %>% gsub('_counts.rds','',.)

# Create directory for output
dir.create(paste0("./Data_generated/", GSE, "/Imputed"), showWarnings = F, recursive = T)

# Iterate over GSM samples and generate rds
for (i in 1:length(GSMs)){
    impute(GSE,GSMs[i])
}

