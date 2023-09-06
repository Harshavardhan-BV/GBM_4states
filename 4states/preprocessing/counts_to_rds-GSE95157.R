source('./common_functions.R')
library(biomaRt)

# Read the counts
df = read.delim("../Data/GSE95157_RAW/GSE95157_counts.txt.gz")
# Convert NCBI gene number to gene symbol
mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesymb = getBM(attributes = c("entrezgene_id","external_gene_name"),filters = 'entrezgene_id', values = df$X, mart = mart)
df = merge(df, genesymb, by.x = "X", by.y = "entrezgene_id")[,-1]
# merge duplicate gene names
df = df %>%
  group_by(external_gene_name) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as gene
rownames(df) = df[,1]
df = df[,-1]
# Remove empty gene name
df = df[!rownames(df) == "",]
rownames(df) = toupper(rownames(df))
# Perform QC and normalization
df = SC_QC(df)
# Save the counts matrix
exportCount(df, 'GSE95157', 'GSE95157', sc=T)
