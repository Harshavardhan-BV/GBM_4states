source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE79671_RAW/GSE79671_CountMatrix.txt.gz", row.names=1)
# geneid to symbol
df = GeneIDToSymb(df)
# TPM to counts
df = countToTpm(df)
# Save the counts matrix
exportCount(df, 'GSE79671', 'GSE79671')
