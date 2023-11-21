source("./common_functions.R")
# Read the counts file
counts = ReadMtx(mtx='../Data/GSE122871_RAW/GSE122871_matrix.mtx.gz', cells='../Data/GSE122871_RAW/GSE122871_barcodes.tsv.gz', features='../Data/GSE122871_RAW/GSE122871_genes.tsv.gz')
# gene name to upper case
rownames(counts) = toupper(rownames(counts))
# Perform QC and normalization
counts = SC_QC(counts)
# Save the counts matrix
exportCount(counts, 'GSE122871', 'GSE122871', sc=T)
