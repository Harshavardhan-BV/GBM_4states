source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE222515_RAW/GSE222515_output_tpm_all_samples.txt.gz",row.names = 1)
# Change gene id to gene names
df = GeneIDToSymb(df)
# Save the counts matrix
exportCount(df, 'GSE222515', 'GSE222515')
