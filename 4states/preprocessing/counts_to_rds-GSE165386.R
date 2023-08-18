source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE165386_RAW/GSE165386_RNAseq_Mouse_raw_counts.txt.gz")
# Set the row names as gene
rownames(df) = df$symbol
# Select only counts
df = df[,1:33]
# to TPM
df = countToTpm(df)
# To uppercase
rownames(df) = toupper(rownames(df))
# Save the counts matrix
exportCount(df, 'GSE165386', 'GSE165386')
