source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE163071_RAW/GSE163071_Raw_Gene_counts.txt.gz")
# Set the row names as gene
rownames(df) = df$Gene.Symbol
df = df[,6:ncol(df)]
# Save the counts matrix
exportCount(df, 'GSE163071', 'GSE163071')
