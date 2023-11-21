source('./common_functions.R')

# List all the files within subdirectories
files <- list.files(path = "../Data/GSE153794_RAW", recursive = TRUE, full.names = TRUE, pattern = "*.txt.gz")

for (file in files){
    print(file)
    # Read the count matrix for all 
    counts = read.delim(file, row.names = 1)
    # merge with the counts matrix
    if (exists("counts_all")){
        counts_all = merge(counts_all,counts[,3,drop=F], by = 0, all = T)
        rownames(counts_all) = counts_all[,1]
        counts_all = counts_all[,-1]
    } else {
        counts_all = counts[,c(2,3)]
    }
}
# merge duplicate genes and take mean of numeric values
counts_all = counts_all %>%
    group_by(gene_short_name) %>%
    summarise_all(mean) %>% data.frame()
# Remove NA genes
counts_all = counts_all[!is.na(counts_all$gene_short_name),]
# Set the row names as geneids
rownames(counts_all) = counts_all[,1]
counts_all = counts_all[,-1]
# Fill na with 0
counts_all[is.na(counts_all)] = 0
# Convert to TPM
counts_all = FPKMToTPM(counts_all)
# Save the counts matrix
exportCount(counts_all, 'GSE153794', 'GSE153794')
