source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE76184_RAW/GSE76184_count_matrix.txt.gz",row.names = 1)
# Convert to TPM
df = countToTpm(df)
# uppercase gene names
rownames(df) = toupper(rownames(df))
# Save the counts matrix
exportCount(df, 'GSE76184', 'GSE76184', sc=T)