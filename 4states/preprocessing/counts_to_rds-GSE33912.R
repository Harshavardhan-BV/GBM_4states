source('./common_functions.R')
library(readxl)
# Read the counts
df = read_excel("../Data/GSE33912_RAW/GSE33912_fpkm.xlsx", sheet=2)
df1 = read_excel("../Data/GSE33912_RAW/GSE33912_fpkm.xlsx", sheet=3, skip = 1)
# Combine the counts
df = cbind(df[,-1],df1[,3:ncol(df1)])
# Uppercase
df$gene_short_name = toupper(df$gene_short_name)
# merge duplicate gene names
df = df %>%
  group_by(gene_short_name) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as gene
rownames(df) = df[,1]
df = df[,-1]
# to TPM
df = FPKMToTPM(df)
# uppercase
# Save the counts matrix
exportCount(df, 'GSE33912', 'GSE33912')
