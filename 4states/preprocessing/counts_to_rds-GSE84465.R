source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE84465_RAW/GSE84465_GBM_All_data.csv.gz",row.names = 1, sep=' ')
# Convert to TPM
df = countToTpm(df)
# Save the counts matrix
exportCount(df, 'GSE84465', 'GSE84465', sc=T)
