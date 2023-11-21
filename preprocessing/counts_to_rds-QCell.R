source('./common_functions.R')

# Read the counts
df = read.delim("../Data/QCell/Normalised-RNA-Sequencing-Data.txt")
# merge duplicate gene names
df = df %>%
  group_by(Gene.Name) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as geneids
rownames(df) = df[,1]
df = df[,-1]
# Convert to TPM
df = FPKMToTPM(df)
# Save the counts matrix
exportCount(df, 'QCell', 'QCell')
