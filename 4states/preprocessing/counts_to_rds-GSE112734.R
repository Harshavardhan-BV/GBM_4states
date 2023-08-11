source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE112734_RAW/GSE112734_AllSamples.GeneExpression.FPKM.txt.gz",row.names = 1)
# Select only counts
df = df[,c(2:32)]
# merge duplicate gene names
df = df %>%
  group_by(Symbol) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as gene
rownames(df) = df[,1]
df = df[,-1]
# Convert to TPM
df = FPKMToTPM(df)
# Save the counts matrix
exportCount(df, 'GSE112734', 'GSE112734')
