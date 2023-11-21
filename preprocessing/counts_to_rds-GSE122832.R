source('./common_functions.R')
library(biomaRt)

# List all the files within subdirectories
files <- list.files(path = "../Data/GSE122832_RAW", recursive = TRUE, full.names = TRUE, pattern = "*.txt.gz")

for (file in files){
    print(file)
    # Read the count matrix for all 
    counts = read.delim(file, row.names=4)
    # Only counts
    counts = counts[,4:ncol(counts)]
    # merge with the counts matrix
    if (exists("counts_all")){
        if (nrow(counts)==nrow(counts_all)){
            counts_all = cbind(counts_all,counts)
        } else {
            print(paste("Error: Row number mismatch", nrow(counts), nrow(counts_all)))
        }
    } else {
        counts_all = counts
    }
}
#
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
genesymb = getBM(attributes = c("refseq_mrna","external_gene_name"),filters = "refseq_mrna", values = rownames(counts_all), mart = mart)
counts_all = merge(counts_all, genesymb, by.x = "row.names", by.y = "refseq_mrna")[,-1]
# merge duplicate gene names
counts_all = counts_all %>%
  group_by(external_gene_name) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as geneids
rownames(counts_all) = counts_all[,1]
counts_all = counts_all[,-1]
# Remove empty gene name
counts_all = counts_all[!rownames(counts_all) == "",]
# To TPM
counts_all = countToTpm(counts_all)
# Save the counts matrix
exportCount(counts_all, 'GSE122832', 'GSE122832')
