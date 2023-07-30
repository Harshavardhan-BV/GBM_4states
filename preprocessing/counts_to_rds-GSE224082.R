source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE224082_RAW/GSE224082_geo_submission_raw_count.txt.gz",row.names = 1)
# geneid to genename
df = GeneIDToSymb(df)
# Convert to TPM
df = countToTpm(df)
# Save the counts matrix
exportCount(df, 'GSE224082', 'GSE224082')