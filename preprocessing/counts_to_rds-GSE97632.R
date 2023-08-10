source('./common_functions.R')
library('readxl')

# Read the counts
files = list.files(path = "../Data/GSE97632_RAW", recursive = TRUE, full.names = TRUE, pattern = "*.xlsx")

for (file in files){
    print(file)
    # Read the count matrix for all 
    counts = read_excel(file)
    # merge with the counts matrix
    if (exists("counts_all")){
        if (nrow(counts)==nrow(counts_all)){
            counts_all = cbind(counts_all, counts[,c(4:ncol(counts))])
        } else {
            print(paste("Error: Row number mismatch", nrow(counts), nrow(counts_all)))
        }
    } else {
        counts_all = counts[,c(2,4:ncol(counts))]

    }
}
# Remove the / from the gene names
counts_all[,1] = gsub("/", "", counts_all[,1])
colnames(counts_all)[1] = "Gene"
# merge duplicate gene names
counts_all = counts_all %>%
  group_by(Gene) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as gene
rownames(counts_all) = counts_all[,1]
counts_all = counts_all[,-1]
# Convert to TPM
counts_all = FPKMToTPM(counts_all)
# Save the counts matrix
exportCount(counts_all, 'GSE97632', 'GSE97632')