source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE98126_RAW/GSE98126_AllSamples.HTSeq.txt.gz",row.names = 1, sep=' ')
# Convert geneid to gene symbol
df = GeneIDToSymb(df)
# Convert to TPM
df = countToTpm(df)
# Save the counts matrix
exportCount(df, 'GSE98126', 'GSE98126')

