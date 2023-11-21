source('./common_functions.R')

# Read the counts
df = read.csv("../Data/GSE179649_RAW/GSE179649_Neuro_Organo_kallisto_gene_tpm.csv.gz",row.names = 1)
# Change gene id to gene names
df = GeneIDToSymb(df)
# Save the counts matrix
exportCount(df, 'GSE179649', 'GSE179649')
