source('./common_functions.R')

# Read the counts
df = read.csv("../Data/GSE154367_RAW/GSE154367_mouseTPM.csv.gz",row.names = 1)
# Convert geneid to symbol
df = GeneIDToSymb(df)
# Convert  to uppercase
rownames(df) = toupper(rownames(df))
# Save the counts matrix
exportCount(df, 'GSE154367', 'GSE154367')