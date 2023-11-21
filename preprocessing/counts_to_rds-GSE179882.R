source('./common_functions.R')

# Read the counts
df = read.csv("../Data/GSE179882_RAW/GSE179882_gexp_counts_sampleTitle.csv.gz", skip=1, row.names=1)
# GeneID to GeneSymbol
df = GeneIDToSymb(df)
# Convert to TPM
df = countToTpm(df)
# Save the counts matrix
exportCount(df, 'GSE179882', 'GSE179882')
