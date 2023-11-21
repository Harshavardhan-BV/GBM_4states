library(dplyr)
library(trqwe)
library(data.table)

df = read.csv("../Data/GSE145556_RAW/GSE145556_RNAseq_mouse_glioblastoma_TPMs.csv.gz")
# Drop rows where df$symbol is NA
df = df[!is.na(df$symbol_human),]
# Select only the counts
df = df[,c(15,2:9)]
# merge duplicate gene names
df = df %>%
  group_by(symbol_human) %>%
  summarize_all(mean) %>% data.frame()
# Set the rownames to the gene names
rownames(df) = df$symbol_human
df = df[,-1]
# Save the data
dir.create("../Data_generated/GSE145556/Counts", showWarnings = FALSE, recursive = TRUE)
mcsaveRDS(df, "../Data_generated/GSE145556/Counts/GSE145556_counts.rds", mc.cores=4)
fwrite(df, "../Data_generated/GSE145556/Counts/GSE145556_counts.tsv", sep="\t", row.names = T)
