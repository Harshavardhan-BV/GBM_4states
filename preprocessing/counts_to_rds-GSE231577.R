source('./common_functions.R')

# Read the counts
df = read.csv("../Data/GSE231577_RAW/GSE231577_tumor_gene_count_matrix.csv.gz",row.names = 1)
# Select only the part after | in colnames
df$genenames = sapply(strsplit(rownames(df), '\\|'), function(x) x[2])
# remove nan genenames
df = df[!is.na(df$genenames),]
# merge duplicate gene names
df = df %>%
  group_by(genenames) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as geneids
rownames(df) = df[,1]
df = df[,-1]
# Convert to TPM
df = countToTpm(df)
# Save the counts matrix
exportCount(df, 'GSE231577-B', 'GSE231577')

# Single cell sample
# Read the ctrl counts
counts_all = Read10X_h5('../Data/GSE231577_RAW/GSM7291138_Ctrl_filtered_feature_bc_matrix.h5')
colnames(counts_all) = paste0('GSM7291138_', colnames(counts_all))
# Read the DIP counts
counts = Read10X_h5('../Data/GSE231577_RAW/GSM7291139_DIP_filtered_feature_bc_matrix.h5')
colnames(counts) = paste0('GSM7291139_', colnames(counts))
# Merge the counts
counts_all = cbind(counts_all, counts)
# Perform QC and normalization
counts_all = SC_QC(counts_all)
# Save the counts matrix
exportCount(counts_all, 'GSE231577', 'GSE231577', sc=T)
