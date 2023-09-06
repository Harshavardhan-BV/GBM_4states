source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE141946_RAW/GSE141946_scRNAseq.counts.txt.gz")
metdat = read.csv("../Data/GSE141946_RAW/GSE141946_scRNAseq.metadata.csv.gz")
# Perform QC and normalization
df = SC_QC(df)
gc(verbose = F, reset=T)
# Save the counts matrix
exportCount(df, 'GSE141946', 'GSE141946', sc=T)
