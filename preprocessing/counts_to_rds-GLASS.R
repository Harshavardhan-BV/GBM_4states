source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GLASS/gene_tpm_matrix_all_samples.tsv",row.names = 1)
# Save the counts matrix
exportCount(df, 'GLASS', 'GLASS')