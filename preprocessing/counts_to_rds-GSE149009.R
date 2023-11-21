source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE149009_RAW/GSE149009_gbm_rna_deseq2.txt.gz",row.names = 1)
# Change gene id to gene names
df = GeneIDToSymb(df)
# Save the counts matrix
exportCount(df, 'GSE149009', 'GSE149009')