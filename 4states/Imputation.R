library(Rmagic)
library(trqwe)

impute <- function(Sample) {
    count_matrix = mcreadRDS(paste0("./Data_generated/",Sample,"_counts.rds"), mc.cores=4)
    count_matrix = t(count_matrix)
    gc()
    count_matrix = magic(count_matrix, n.jobs = -2, solver='exact')
    gc()
    mcsaveRDS(t(count_matrix$result), paste0("./Data_generated/",Sample,"_imputed.rds"), mc.cores = 4)
}

# GSE ID
GSE = "GSE168004"

# GSM file list
GSMs = c("OSM_celllines", "mgg23", "mgg75")

# Iterate over GSM samples and generate rds
for (i in 1:length(GSMs)){
    impute(paste0(GSE,'/',GSMs[i]))
}

