library(dplyr)
library(trqwe)
library(data.table)
# List all the files within subdirectories of ./Data/TCGA
files <- list.files(path = "../Data/TCGA", recursive = TRUE, full.names = TRUE, pattern = "*.tsv")

for (file in files){
    print(file)
    # Read the count matrix for all 
    counts = read.delim(file,skip=1)
    # Remove the first 5 rows and select only the gene_name and tpm_unstranded columns
    counts = counts[5:nrow(counts),c("gene_name","tpm_unstranded")]
    # merge with the counts matrix
    if (exists("counts_all")){
        if (nrow(counts)==nrow(counts_all)){
            counts_all = cbind(counts_all,counts$tpm_unstranded)
        } else {
            print(paste("Error: Row number mismatch", nrow(counts), nrow(counts_all)))
        }
    } else {
        counts_all = counts
    }
}
# Set the column names
colnames(counts_all)[2:ncol(counts_all)] = seq(1,ncol(counts_all)-1)
# Merge duplicate row names and take mean
counts_all = counts_all %>%
  group_by(gene_name) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as gene names
rownames(counts_all) = counts_all$gene_name
counts_all = counts_all[,-1]
# Export the counts for the GBM celllines only as tsv
dir.create('../Data_generated/TCGA/Counts', showWarnings = FALSE, recursive = TRUE)
fwrite(counts_all, '../Data_generated/TCGA/Counts/TCGA_counts.tsv', sep='\t', row.names = T)
# Export the counts for GBM celllines as rds
mcsaveRDS(counts_all, "../Data_generated/TCGA/Counts/TCGA_counts.rds", mc.cores=4)
