source('./common_functions.R')

# List all the files within subdirectories
files <- list.files(path = "../Data/GSE147352_RAW", recursive = TRUE, full.names = TRUE, pattern = "*.txt.gz")

for (file in files){
    print(file)
    # Read the count matrix for all 
    counts = read.delim(file, header=F, row.names=1)
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
# Select the Sample names from the file names
colnames(counts_all) = gsub("_.*", "", basename(files))
# GeneID to GeneName
counts_all = GeneIDToSymb(counts_all)
# Convert to TPM
counts_all = countToTpm(counts_all)
# Save the counts matrix
exportCount(counts_all, 'GSE147352', 'GSE147352')
