source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE218042_RAW/GSE218042_CCR7shRNA_TPM.txt.gz")
# merge duplicate gene names
df = df %>%
  group_by(Gene) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as gene
rownames(df) = df[,1]
df = df[,-1]
# Save the counts matrix
exportCount(df, 'GSE218042', 'GSE218042')
