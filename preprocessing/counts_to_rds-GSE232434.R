source('./common_functions.R')
library(biomaRt)

# Read the count matrix
df = read.csv("../Data/GSE232434_RAW/GSE232434_MIFsnp.hg19KnownGene.tpm.csv.gz", row.names = 1)
# Convert NCBI gene number to gene symbol
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genesymb = getBM(attributes = c("entrezgene_id","external_gene_name"),filters = 'entrezgene_id', values = rownames(df), mart = mart)
df = merge(df, genesymb, by.x = "row.names", by.y = "entrezgene_id")[,-1]
# merge duplicate gene names
df = df %>%
  group_by(external_gene_name) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as geneids
rownames(df) = df[,1]
df = df[,-1]
# Remove empty gene name
df = df[!rownames(df) == "",]
# Save the counts matrix
exportCount(df, 'GSE232434', 'GSE232434')