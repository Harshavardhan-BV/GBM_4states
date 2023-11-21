source('./common_functions.R')
library(readxl)

# Read the counts
df = read_excel("../Data/GSE117632_RAW/GSE117632_Hypoxia_ALKBH1_Knockdown_RNAseq_Processed.xlsx")
# Convert to dataframe
df = data.frame(df)
# Set the row names as gene
rownames(df) = df[,1]
df = df[,-1]
# Save the counts matrix
exportCount(df, 'GSE117632', 'GSE117632')
