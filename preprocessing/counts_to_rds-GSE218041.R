source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE218041_RAW/GSE218041_CCR7_inhibitor_TPM.txt.gz")
# merge duplicate gene names
df = df %>%
  group_by(Gene) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as gene
rownames(df) = df[,1]
df = df[,-1]
# Save the counts matrix
exportCount(df, 'GSE218041', 'GSE218041')
