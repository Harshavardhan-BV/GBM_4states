source("./common_functions.R")
# Read the seurat object
counts = Read10X_h5('../Data/GSE213308_RAW/GSE213308_filtered_feature_bc_matrix.h5')
# gene name to upper case
rownames(counts) = toupper(rownames(counts))
# Perform QC and normalization
counts = SC_QC(counts)
# Save the counts matrix
exportCount(counts, 'GSE213308', 'GSE213308', sc=T)
