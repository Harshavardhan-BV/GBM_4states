source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE200062_RAW/GSE200062_Human_Normalized_counts.tsv.gz")
# Remove na genes
df = df[!is.na(df$gene_name),-1]
# merge duplicate gene names
df = df %>%
  group_by(gene_name) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as gene
rownames(df) = df[,1]
df = df[,-1]
# Save the counts matrix
exportCount(df, 'GSE200062', 'Human')


# # Read the counts
df = read.delim("../Data/GSE200062_RAW/GSE200062_Mouse_Normalized_counts.tsv.gz")
# Remove na genes
df = df[!is.na(df$external_gene_name),-1]
# to uppercase
df$external_gene_name = toupper(df$external_gene_name)
# merge duplicate gene names
df = df %>%
  group_by(external_gene_name) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as gene
rownames(df) = df[,1]
df = df[,-1]
# Save the counts matrix
exportCount(df, 'GSE200062', 'Mouse')