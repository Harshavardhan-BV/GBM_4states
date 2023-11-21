source('./common_functions.R')
library(readxl)

# Read the counts
df = read_excel('../Data/GSE135408_RAW/GSE135408_Tissue_raw_counts.xlsx')
# merge duplicate gene names
df = df %>% group_by(EnsemblID) %>% summarise_all(sum) %>% data.frame()
# set rownames
rownames(df) = df$EnsemblID
df = df[,-1]
# Convert rownames to gene symbols
df = GeneIDToSymb(df)
# to TPM
df = countToTpm(df)
# Save the counts matrix
exportCount(df, 'GSE135408', 'Tissue')


# Read the counts
df = read_excel('../Data/GSE135408_RAW/GSE135408_Clone_raw_counts.xlsx')
# merge duplicate gene names
df = df %>% group_by(EnsemblID) %>% summarise_all(sum) %>% data.frame()
# set rownames
rownames(df) = df$EnsemblID
df = df[,-1]
# Convert rownames to gene symbols
df = GeneIDToSymb(df)
# to TPM
df = countToTpm(df)
# Save the counts matrix
exportCount(df, 'GSE135408', 'Clone')

# Read the counts
df = read_excel('../Data/GSE135408_RAW/GSE135408_CD38_raw_counts.xlsx')
# merge duplicate gene names
df = df %>% group_by(EnsemblID) %>% summarise_all(sum) %>% data.frame()
# set rownames
rownames(df) = df$EnsemblID
df = df[,-1]
# Convert rownames to gene symbols
df = GeneIDToSymb(df)
# to TPM
df = countToTpm(df)
# uppercase
rownames(df) = toupper(rownames(df))
# Save the counts matrix
exportCount(df, 'GSE135408', 'CD38')
