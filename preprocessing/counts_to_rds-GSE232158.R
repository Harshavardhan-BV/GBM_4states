source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE232158_RAW/GSE232158_mRNA_Expression_submit.txt.gz",row.names = 1)
# additional samples?
df1 = read.delim("../Data/GSE232158_RAW/GSE232158_mRNA_Expression_submit2.tsv.gz",row.names = 1)
# Select only the counts
df = df[,c(6:ncol(df))]
df1 = df1[,c(6:ncol(df1))]
# merge the two dataframes
df = cbind(df,df1)
# geneid to gene_symbol
df = GeneIDToSymb(df)
# Convert to TPM
df = countToTpm(df)
# Save the counts matrix
exportCount(df, 'GSE232158', 'GSE232158')

