library(dplyr)
library(trqwe)

df = read.delim("../Data/GSE57872_RAW/GSE57872_GBM_data_matrix.txt.gz",row.names = 1)
dir.create("../Data_generated/GSE57872/Counts", showWarnings = FALSE, recursive = TRUE)
mcsaveRDS(df, "../Data_generated/GSE57872/Counts/GSE57872_counts.rds", mc.cores=4)
