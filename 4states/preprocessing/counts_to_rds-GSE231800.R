source("./common_functions.R")
library('Seurat')
# Read the seurat object
counts = mcreadRDS('../Data/GSE231800_RAW/GSE231800_S4F_50.rds', mc.cores=4)
# Extract counts 
counts = GetAssayData(object = counts, assay='RNA')
gc(verbose = F, reset=T)
# gene name to upper case
rownames(counts) = toupper(rownames(counts))
# Perform QC and normalization
counts = SC_QC(counts)
# Save the counts matrix
exportCount(counts, 'GSE231800', 'GSE231800', sc=T)
