source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE139533_RAW/GSE139533_mergede_protein_coding_counts.tab.gz",row.names = 1)
# Convert to TPM
df = countToTpm(df)
# Save the counts matrix
exportCount(df, 'GSE139533', 'GSE139533')
