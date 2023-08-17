source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE138662_RAW/GSE138662_T98G_time_course_gene_counts.tsv.gz", row.names=1)
# gene id to symbol
df = GeneIDToSymb(df)
# to TPM
df = countToTpm(df)
# Save the counts matrix
exportCount(df, 'GSE138662', 'GSE138662')
