source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE184941_RAW/GSE184941_Raw_gene_counts_matrix.txt.gz",row.names = 1)
# Change gene id to gene names
df = GeneIDToSymb(df)
# Convert to tpm
df = countToTpm(df)
# Save the counts matrix
exportCount(df, 'GSE184941', 'GSE184941')
