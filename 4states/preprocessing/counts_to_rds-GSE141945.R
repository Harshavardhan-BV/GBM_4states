source('./common_functions.R')

# Read the counts
df = read.csv("../Data/GSE141945_RAW/GSE141945_RNAseq.counts.csv.gz")
# merge duplicate gene names
df = df %>%
  group_by(X) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as geneids
rownames(df) = df[,1]
df = df[,-1]
# Convert to TPM
df = countToTpm(df)
# Save the counts matrix
exportCount(df, 'GSE141945', 'GSE141945',)
