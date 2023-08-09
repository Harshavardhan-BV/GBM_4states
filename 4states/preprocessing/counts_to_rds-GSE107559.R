source('./common_functions.R')
# Read the counts
df = read.csv("../Data/GSE107559_RAW/GSE107559_ivygap_fpkm_table.csv.gz",row.names = 1)
# Read the gene names
genes = read.csv("../Data/GSE107559_RAW/GSE107559_ivygap_rows-genes.csv.gz", row.names=1)
# Replace the row names with mapping from genes
rownames(df) = genes[rownames(df),"gene_symbol"]
# RPKM to TPM
df = FPKMToTPM(df)
# Save the counts matrix
exportCount(df, 'GSE107559', 'GSE107559')
