source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE213307_RAW/GSE213307_RNASeq.Counts.tsv.gz",row.names = 1)
# Convert to TPM
df = countToTpm(df)
# uppercase gene names
rownames(df) = toupper(rownames(df))
# Save the counts matrix
exportCount(df, 'GSE213307', 'GSE213307')