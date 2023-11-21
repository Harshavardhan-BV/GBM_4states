source('./common_functions.R')

# Read the data
df = read.csv("../Data/GSE212246_RAW/GSE212246_TPM_qc.csv.gz",row.names=1)
metdat = read.csv("../Data/GSE212246_RAW/GSE212246_metadata_mGBM.csv.gz",row.names=1)
# Convert rownames to uppercase
rownames(df) = toupper(rownames(df))
# Select only the malignant samples
gcells = rownames(metdat[metdat$annotate=="Malignant",])
df = df[,gcells]
# QC
df = SC_QC(df, normalize = F)
# Save the data
exportCount(df, 'GSE212246','GSE212246',sc=T)
