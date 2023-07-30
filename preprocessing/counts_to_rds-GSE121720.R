source('./common_functions.R')

# Read the count matrix
df = read.delim("../Data/GSE121720_RAW/GSE121720_RNAseq_expression_matrix_TPMs.txt.gz", row.names = 1)

# Convert gene IDs to gene symbols
df = GeneIDToSymb(df)

# Save the counts matrix
exportCount(df, 'GSE121720', 'GSE121720')