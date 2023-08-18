source('./common_functions.R')

# List all the files within subdirectories
files <- list.files(path = "../Data/GSE155434_RAW", recursive = TRUE, full.names = TRUE, pattern = "*_feature_count.txt.gz")

for (file in files){
    print(file)
    # Read the count matrix for all 
    counts = read.delim(file, skip=1, row.names=1)
    # Only counts
    counts = counts[,6,drop=F]
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
# Change geneid to genename
counts_all = GeneIDToSymb(counts_all)
# to TPM
counts_all = countToTpm(counts_all)
# Save the counts matrix
exportCount(counts_all, 'GSE155434', 'GSE155434')
