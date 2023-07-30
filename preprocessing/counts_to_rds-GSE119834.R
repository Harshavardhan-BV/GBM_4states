source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE119834_RAW/GSE119834_fpkm_table.txt.gz",row.names = 1)
# Save the counts matrix
exportCount(df, 'GSE119834', 'GSE119834')
