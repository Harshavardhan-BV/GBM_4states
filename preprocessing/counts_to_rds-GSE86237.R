source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE86237_RAW/GSE86237_geo_matrix.txt.gz")
# merge duplicate gene names
df = df[,-2] %>%
  group_by(gene_id) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as gene
rownames(df) = df[,1]
df = df[,-1]
# Convert to TPM
df = FPKMToTPM(df)
# Save the counts matrix
exportCount(df, 'GSE86237', 'GSE86237')
