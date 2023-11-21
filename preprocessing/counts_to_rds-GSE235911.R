source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE235911_RAW/GSE235911_bulkRNAseq_count_table.txt.gz",row.names = 1)
# select only the counts columns
df = df[,c(6:ncol(df))]
# Convert to TPM
df = countToTpm(df)
# uppercase gene names
rownames(df) = toupper(rownames(df))
# Save the counts matrix
exportCount(df, 'GSE235911', 'GSE235911')
