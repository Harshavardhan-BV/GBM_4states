source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE114222_RAW/GSE114222_geneExpressionTable.txt.gz", row.names=1)
# Select only counts
df = df[,-ncol(df)]
# Convert to TPM
df = countToTpm(df)
# Save the counts matrix
exportCount(df, 'GSE114222', 'GSE114222')
