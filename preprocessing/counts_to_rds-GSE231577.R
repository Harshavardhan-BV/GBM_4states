source('./common_functions.R')

# Read the counts
df = read.csv("../Data/GSE231577_RAW/GSE231577_tumor_gene_count_matrix.csv.gz",row.names = 1)
# Select only the part after | in colnames
df$genenames = sapply(strsplit(rownames(df), '\\|'), function(x) x[2])
# remove nan genenames
df = df[!is.na(df$Name),]
# merge duplicate gene names
df = df %>%
  group_by(genenames) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as geneids
rownames(df) = df[,1]
df = df[,-1]
# Convert to TPM
df = countToTpm(df)
# Save the counts matrix
exportCount(df, 'GSE231577-B', 'GSE231577')
