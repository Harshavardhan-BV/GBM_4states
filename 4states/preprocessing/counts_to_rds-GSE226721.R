source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE226721_RAW/GSE226721_gene_count.txt.gz",row.names = 1)
df = df[,c(1:9)]
# merge duplicate rows
df = df %>%
  group_by(gene_name) %>%
  summarize_all(mean) %>% data.frame()
# set rownames
rownames(df) = df$gene_name
df = df[,-1]
# Save the counts matrix
exportCount(df, 'GSE226721', 'GSE226721')
