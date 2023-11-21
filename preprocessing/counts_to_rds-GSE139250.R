source('./common_functions.R')

# List all the files within subdirectories
files <- list.files(path = "../Data/GSE139250_RAW", recursive = TRUE, full.names = TRUE, pattern = ".RCC.gz")

for (file in files){
    print(file)
    # Read the count matrix for all 
    counts = read.csv(file, skip = 26)
    # Select the columns with geneid and counts
    counts = counts[,c(2,4)]
    # merge with the counts matrix
    if (exists("counts_all")){
        if (nrow(counts)==nrow(counts_all)){
            counts_all = cbind(counts_all,counts[,2])
        } else {
            print(paste("Error: Row number mismatch", nrow(counts), nrow(counts_all)))
        }
    } else {
        counts_all = counts
    }
}
# Select the part b4 _ as sample names
colnames(counts_all)[-1] = gsub('_.*','',basename(files))
# remove nan genenames
counts_all = counts_all[!is.na(counts_all$Name),]
# merge duplicate gene names
counts_all = counts_all %>%
  group_by(Name) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as geneids
rownames(counts_all) = counts_all[,1]
counts_all = counts_all[,-1]

# Convert to TPM
counts_all = countToTpm(counts_all)
# Save the counts matrix
exportCount(counts_all, 'GSE139250', 'GSE139250')
