library(dplyr)
library(trqwe)

# Read the data
df = read.csv("../Data/GSE212246_RAW/GSE212246_TPM_qc.csv.gz",row.names=1)
metdat = read.csv("../Data/GSE212246_RAW/GSE212246_metadata_mGBM.csv.gz",row.names=1)
# Convert rownames to uppercase
rownames(df) = toupper(rownames(df))
# Select only the malignant samples
gcells = rownames(metdat[metdat$annotate=="Malignant",])
df = df[,gcells]
# Save the data
dir.create("../Data_generated/GSE212246/Counts", showWarnings = FALSE, recursive = TRUE)
mcsaveRDS(df, "../Data_generated/GSE212246/Counts/GSE212246_counts.rds", mc.cores=4)
