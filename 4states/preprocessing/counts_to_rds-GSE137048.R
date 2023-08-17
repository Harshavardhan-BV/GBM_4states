source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE137048_RAW/GSE137048_counts.txt.gz", row.names=1)
# Gene ID to symbol
df = GeneIDToSymb(df)
# Convert to TPM
df = countToTpm(df)
# Save the counts matrix
exportCount(df, 'GSE137048', 'GSE137048')
