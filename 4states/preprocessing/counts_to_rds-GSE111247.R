source('./common_functions.R')

# Read the counts
df = read.csv("../Data/GSE111247_RAW/GSE111247_all_counts.csv.gz", row.names=1)
# Convert to TPM
df = countToTpm(df)
# Save the counts matrix
exportCount(df, 'GSE111247', 'GSE111247')
