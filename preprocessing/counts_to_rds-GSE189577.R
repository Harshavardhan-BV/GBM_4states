source('./common_functions.R')

# Read the counts
counts_all = read.csv("../Data/GSE189577_RAW/GSM5704141_EarlyStage_Day1.csv.gz",row.names = 1)
counts = read.csv("../Data/GSE189577_RAW/GSM5704142_EarlyStage_Day2.csv.gz",row.names = 1)
# Merge the counts
counts_all = merge(counts_all, counts, by = 'row.names', all = T)
# Set rownames
rownames(counts_all) = counts_all[,1]
counts_all = counts_all[,-1]
# Set nan to 0
counts_all[is.na(counts_all)] = 0
# Convert to TPM
counts_all = SC_QC(counts_all)
# uppercase gene names
rownames(counts_all) = toupper(rownames(counts_all))
# Save the counts matrix
exportCount(counts_all, 'GSE189577', 'GSE189577', sc=T)