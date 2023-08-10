source('./common_functions.R')
library('readxl')

# Read the counts
df = read_excel("../Data/GSE154958_RAW/GSE154958_cell_culture_salmon_tpm.xlsx")
df = data.frame(df[,-1])
# merge duplicate gene names
df = df %>%
  group_by(Gene.Symbol) %>%
  summarize_all(mean) %>% data.frame()
df = df[!is.na(df$Gene.Symbol),]
# Set rownames
rownames(df) = df[,1]
df = df[,-1]
# Save the counts matrix
exportCount(df, 'GSE154958', 'cellculture')

# Read the counts
df = read_excel("../Data/GSE154958_RAW/GSE154958_ffpe_salmon_tpm_all_samples.xlsx")
df = data.frame(df[,-1])
# merge duplicate gene names
df = df %>%
  group_by(Gene.Symbol) %>%
  summarize_all(mean) %>% data.frame()
df = df[!is.na(df$Gene.Symbol),]
# Set rownames
rownames(df) = df[,1]
df = df[,-1]
# Save the counts matrix
exportCount(df, 'GSE154958', 'ffpe')