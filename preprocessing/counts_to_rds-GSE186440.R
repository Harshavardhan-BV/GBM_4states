source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE186440_RAW/GSE186440_DAXX_ATRX_RNAseq.txt.gz", row.names=1)
# GeneID to GeneSymbol
df = GeneIDToSymb(df)
# Convert to TPM
df = countToTpm(df)
# Save the counts matrix
exportCount(df, 'GSE186440', 'GSE186440')