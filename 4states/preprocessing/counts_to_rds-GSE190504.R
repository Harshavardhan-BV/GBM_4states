source('./common_functions.R')
library('readxl')

# Read the counts
df = read_excel("../Data/GSE190504_RAW/GSE190504_Processed_Data_Spreadsheet_Glioma_Study.xlsx", skip = 4)
# Transpose
df = t(df)
# Set the column name
colnames(df) = df[1,]
df = df[-c(1,2),]
# Store the rowname
rownam = rownames(df)
# Convert to numeric dataframe
df = data.frame(apply(df, 2, as.numeric))
rownames(df) = rownam
# Convert to TPM
df = CPMToFPKM(df)
df = FPKMToTPM(df)
# Save the counts matrix
exportCount(df, 'GSE190504', 'GSE190504')

