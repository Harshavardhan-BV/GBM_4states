source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE77530_RAW/GSE77530_GBM_AH_32_RSEQ_expression_profile.txt.gz", row.names=1)
# To TPM
df = FPKMToTPM(df)
# Save the counts matrix
exportCount(df, 'GSE77530', 'GSE77530')
