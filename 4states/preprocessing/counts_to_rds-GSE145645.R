source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE145645_RAW/GSE145645_XLGBM_FPKMlg2.txt.gz",row.names = 1)
df = df[,-2]
# merge duplicate gene names
df = df %>%
  group_by(Gene) %>%
  summarize_all(mean) %>% data.frame()
# set rownames
rownames(df) = df[,1]
df = df[,-1]
# Convert to TPM
df = FPKMToTPM(df)
# Save the counts matrix
exportCount(df, 'GSE145645', 'GSE145645')
