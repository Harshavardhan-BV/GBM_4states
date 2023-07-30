library(dplyr)
library(trqwe)
library(data.table)
source("common_functions.R")

df = read.delim("../Data/GSE176539_RAW/GSE176539_GSC_TS-543_DMSO_HLM006474_S3I-201_combin_TPM.txt.gz",row.names = 1)
# Select only the columns with counts
df = df[,c(1,3:10)]
# remove the nan genes
df = df[!(df$name=='nan'),]
# 
df = df %>%
  group_by(name) %>%
  summarize_all(mean) %>% data.frame()
# set row names as gene names
rownames(df) = df$name
df = df[,-1]
# Save the data
dir.create("../Data_generated/GSE176539/Counts", showWarnings = FALSE, recursive = TRUE)
mcsaveRDS(df, "../Data_generated/GSE176539/Counts/GSE176539_counts.rds", mc.cores=4)
fwrite(df, "../Data_generated/GSE176539/Counts/GSE176539_counts.tsv", sep="\t", row.names = T)
