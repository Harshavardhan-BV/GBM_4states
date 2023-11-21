source('./common_functions.R')

# List all the files within subdirectories
files <- list.files(path = "../Data/GSE183910_RAW", recursive = TRUE, full.names = TRUE, pattern = "*.txt.gz")

for (file in files){
    print(file)
    # Read the count matrix for all 
    counts = read.delim(file)
    # merge with the counts matrix
    if (exists("counts_all")){
        if (nrow(counts)==nrow(counts_all)){
            counts_all = cbind(counts_all,counts[,-1])
        } else {
            print(paste("Error: Row number mismatch", nrow(counts), nrow(counts_all)))
        }
    } else {
        counts_all = counts
    }
}
# Set the row names as geneids
# merge duplicate gene names
counts_all = counts_all %>%
  group_by(X) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as gene
rownames(counts_all) = counts_all[,1]
counts_all = counts_all[,-1]
# Convert to TPM
counts_all = CPMToFPKM(counts_all)
counts_all = FPKMToTPM(counts_all)
# Convert to uppercase
rownames(counts_all) = toupper(rownames(counts_all))
# Save the counts matrix
exportCount(counts_all, 'GSE183910', 'GSE183910')
