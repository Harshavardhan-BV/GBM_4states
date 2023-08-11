source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE156819_RAW/GSE156819_m6_f6_dmso_cpi_RNAseq_all_count.txt.gz",row.names = 1)
# Convert geneid to gene name
df = GeneIDToSymb(df)
# Convert to TPM
df = countToTpm(df)
# uppercase gene names
rownames(df) = toupper(rownames(df))
# Save the counts matrix
exportCount(df, 'GSE156819', 'cpi')

# Read the counts
df = read.delim("../Data/GSE156819_RAW/GSE156819_m6_f6_dmso_jq1_RNAseq_all_count_sortID.txt.gz",row.names = 1)
# merge duplicate gene names
df = df %>%
  group_by(gene_symbol) %>%
  summarize_all(mean) %>% data.frame()
# drop na genes
df = df[!is.na(df$gene_symbol),]
# Set the row names as geneids
rownames(df) = df[,1]
df = df[,-1]
# Convert to TPM
df = countToTpm(df)
# uppercase gene names
rownames(df) = toupper(rownames(df))
# Save the counts matrix
exportCount(df, 'GSE156819', 'jq1')