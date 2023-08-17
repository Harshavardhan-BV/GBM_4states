source('./common_functions.R')

# List all the files within subdirectories
files <- list.files(path = "../Data/GSE198438_RAW", recursive = TRUE, full.names = TRUE, pattern = "*.sf.gz")

for (file in files){
    print(file)
    # Read the count matrix for all 
    counts = read.delim(file, row.names=1)
    # Only counts
    counts = counts[,3,drop=F]
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
# Set the column names as sample names
colnames(counts_all) = gsub('_.*','',basename(files))
# Change geneid to genename
counts_all = GeneIDToSymb(counts_all)
# uppercase
rownames(counts_all) = toupper(rownames(counts_all))
# Save the counts matrix
exportCount(counts_all, 'GSE198438', 'GSE198438')
