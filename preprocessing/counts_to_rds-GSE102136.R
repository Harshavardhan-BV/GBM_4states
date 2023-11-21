source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE102136_RAW/GSE102136_Read_counts_by_clone_EnsemblID.txt.gz",row.names = 1)
# geneid to gene_symbol
df = GeneIDToSymb(df)
# Convert to TPM
df = countToTpm(df)
# Save the counts matrix
exportCount(df, 'GSE102136', 'GSE102136')
