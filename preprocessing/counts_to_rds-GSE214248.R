source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE214248_RAW/GSE214248_core_table_gene.txt.gz",row.names = 1)
# Select only the tpm
df = df[,c(2,7:10)]
# Set row names as gene names
rownames(df) = df[,1]
df = df[,-1]
# Save the counts matrix
exportCount(df, 'GSE189577', 'GSE214248')