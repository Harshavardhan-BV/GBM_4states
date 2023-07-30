source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE154041_RAW/GSE154041_samples_expression_table_filtered.txt.gz",row.names = 1)
# drop na rows
df = df[!is.na(df$Gene_name),]
# Select all but the last 2 columns
df = df[,-c(ncol(df),ncol(df)-1)]
# merge duplicate gene names
df = df %>%
  group_by(Gene_name) %>%
  summarize_all(mean) %>% data.frame()
# Set the rownames to the gene names
rownames(df) = df$Gene_name
df = df[,-1]
# Save the counts matrix
exportCount(df, 'GSE154041', 'GSE154041')
