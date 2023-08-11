source('./common_functions.R')

# List all the files within subdirectories
files <- list.files(path = "../Data/GSE181144_RAW", recursive = TRUE, full.names = TRUE, pattern = "*.txt.gz")

for (file in files){
    print(file)
    # Read the count matrix for all 
    counts = read.delim(file, row.names=1)
    # Reverse stranded
    counts = counts[,6, drop=F]
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
# GeneID to GeneName
counts_all = GeneIDToSymb(counts_all)
# Convert to TPM
counts_all = countToTpm(counts_all)
# Save the counts matrix
exportCount(counts_all, 'GSE181144', 'GSE181144')
