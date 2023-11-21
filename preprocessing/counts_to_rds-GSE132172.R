source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE132172_RAW/GSE132172_glio_chrom_instability_normalized_expression_matrix.tsv.gz")
# Perform QC and normalization
df = SC_QC(df, normalize = F)
# Save the counts matrix
exportCount(df, 'GSE132172', 'GSE132172', sc=T)
