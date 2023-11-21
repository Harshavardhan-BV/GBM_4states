source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE158550_RAW/GSE158550_All_FPKM_normalized_combined.txt.gz",row.names = 1)
# Convert to TPM
df = FPKMToTPM(df)
# Remove df$Xs. WTF is it even????
df = df[,!grepl("X",colnames(df))]
# Save the counts matrix
exportCount(df, 'GSE158550', 'GSE158550')
