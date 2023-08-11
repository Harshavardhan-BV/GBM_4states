source('./common_functions.R')
library(readxl)

# Read the counts
df = read_excel("../Data/GSE112735_RAW/GSE112735_mRNA_Expression_Profiling.xlsx", skip=18)
# Convert to dataframe
df = data.frame(df[,c(1,7:12)])
# Merge duplicate gene names
df = df %>%
  group_by(gene_short_name) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as gene
rownames(df) = df[,1]
df = df[,-1]
# Convert to TPM
df = FPKMToTPM(df)
# Save the counts matrix
exportCount(df, 'GSE112735', 'GSE112735')
