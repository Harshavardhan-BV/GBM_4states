source('./common_functions.R')

# Read counts
counts = read.delim("../Data/GSE57872_RAW/GSE57872_GBM_data_matrix.txt.gz",row.names = 1)
# QC
MT = grep(pattern = "^MT-", x = rownames(counts), value = TRUE)
percent.mt = colSums(counts[MT, ])/colSums(counts)*100
# Subset for outliers using 5MADs for total counts and 3MADs for percent.mt
counts = counts[, is_not_outlier(total.counts, 5) & is_not_outlier(total.counts, -5)  & is_not_outlier(percent.mt, 3)]
# Save counts
exportCount(counts, "GSE57872", "GSE57872", sc=T)
