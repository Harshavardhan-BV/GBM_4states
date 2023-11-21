source('./common_functions.R')

# Read the counts
df = read.csv("../Data/GSE189517_RAW/GSE189517_mouseGliomasPDGFB_expected_counts.csv.gz",row.names = 1)
# Remove na from gene names
df = df[!is.na(df$external_gene_name),]
# merge duplicate gene names
df = df %>%
  group_by(external_gene_name) %>%
  summarize_all(mean) %>% data.frame()
# Set the rownames to the gene names
rownames(df) = df$external_gene_name
df = df[,-1]
# Convert to TPM
df = countToTpm(df)
# uppercase gene names
rownames(df) = toupper(rownames(df))
# Save the counts matrix
exportCount(df, 'GSE189517', 'GSE189517')

