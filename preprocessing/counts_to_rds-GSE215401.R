source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE215401_RAW/GSE215401_GSC11_V_CM.txt.gz", row.names=3)
# merge duplicate gene names
df = df[3:ncol(df)]
# Convert to TPM
df = FPKMToTPM(df)
# Save the counts matrix
exportCount(df, 'GSE215401', 'V_CM')

# Read the counts
df = read.delim("../Data/GSE215401_RAW/GSE215401_V_C_M_CM.txt.gz")
# Merge duplicate gene names
df = df %>% group_by(Symbol) %>% summarize_all(mean) %>% data.frame()
# Set the row names as gene
rownames(df) = df[,1]
df = df[,-1]
# Convert to TPM
df = FPKMToTPM(df)
# Save the counts matrix
exportCount(df, 'GSE215401', 'V_C_M_CM')
