library(dplyr)
library(trqwe)

df = read.delim("./Data/GSE131928_RAW/GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz",row.names = 1)
mcsaveRDS(df, "./Data_generated/GSE131928/GSM3828672_counts.rds", mc.cores=4)

df = read.delim("./Data/GSE131928_RAW/GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv.gz",row.names = 1)
mcsaveRDS(df, "./Data_generated/GSE131928/GSM3828673_counts.rds", mc.cores=4)
