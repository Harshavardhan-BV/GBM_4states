source('./common_functions.R')

df = read.delim("../Data/GSE131928_RAW/GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz",row.names = 1)
df = SC_QC(df,normalize = F)
exportCount(df, "GSE131928", "GSM3828672", sc=T)

df = read.delim("../Data/GSE131928_RAW/GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv.gz",row.names = 1)
df = SC_QC(df,normalize = F)
exportCount(df, "GSE131928", "GSM3828673", sc=T)
