source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE178337_RAW/GSE178337_CCL21_TPM.txt.gz",row.names = 1)
# Save the counts matrix
exportCount(df, 'GSE178337', 'GSE178337')