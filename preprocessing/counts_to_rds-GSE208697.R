source('./common_functions.R')

# List all the files within subdirectories
files <- list.files(path = "../Data/GSE208697_RAW", recursive = TRUE, full.names = TRUE, pattern = "*.txt.gz")

for (file in files){
    print(file)
    # Read the count matrix for all 
    counts = read.delim(file,row.names = 1)
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

# Gene ID to Gene Symbol
counts_all = GeneIDToSymb(counts_all)
# convert to TPM
counts_all = countToTpm(counts_all)
# save the counts matrix
exportCount(counts_all, 'GSE208697', 'GSE208697', sc=FALSE)