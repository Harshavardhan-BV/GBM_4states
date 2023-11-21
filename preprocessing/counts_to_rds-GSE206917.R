source('./common_functions.R')
library(readxl)
# Read the counts
df = read_excel("../Data/GSE206917_RAW/GSE206917_tissue_gene_expression.xls",sheet=1)
# Select only the columns with counts
df = df[,c(3:ncol(df))]
# merge duplicate rows
df = df %>%
  group_by(gene_symbol) %>%
  summarize_all(mean) %>% data.frame()
# set rownames
rownames(df) = df$gene_symbol
df = df[,-1]
# Save the counts matrix
exportCount(df, 'GSE206917', 'tissue')

# Read the counts
df = read_excel("../Data/GSE206917_RAW/GSE206917_U87_gene_expression.xls",sheet=1)
# Select only the columns with counts
df = df[,c(3:ncol(df))]
# merge duplicate rows
df = df %>%
  group_by(gene_symbol) %>%
  summarize_all(mean) %>% data.frame()
# set rownames
rownames(df) = df$gene_symbol
df = df[,-1]
# Save the counts matrix
exportCount(df, 'GSE206917', 'U87')

# Read the counts
df = read_excel("../Data/GSE206917_RAW/GSE206917_U251_gene_expression.xls",sheet=1)
# Select only the columns with counts
df = df[,c(3:ncol(df))]
# merge duplicate rows
df = df %>%
  group_by(gene_symbol) %>%
  summarize_all(mean) %>% data.frame()
# set rownames
rownames(df) = df$gene_symbol
df = df[,-1]
# Save the counts matrix
exportCount(df, 'GSE206917', 'U251')

