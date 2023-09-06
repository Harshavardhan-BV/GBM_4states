source('./common_functions.R')

# Read the counts
df = read.csv("../Data/GSE195682_RAW/GSE195682_all_count_matrix.csv.gz",row.names = 1)
# QC and norm
df = SC_QC(df)
# Save the counts matrix
exportCount(df, 'GSE195682', 'GSE195682', sc=T)
