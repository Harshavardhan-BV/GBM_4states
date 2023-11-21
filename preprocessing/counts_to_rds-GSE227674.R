source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE227674_RAW/GSE227674_raw_counts.txt.gz",row.names = 1)
# Convert geneid to genename
df = GeneIDToSymb(df)
# Convert to TPM
df = countToTpm(df)
# Save the counts matrix
exportCount(df, 'GSE227674', 'GSE227674')