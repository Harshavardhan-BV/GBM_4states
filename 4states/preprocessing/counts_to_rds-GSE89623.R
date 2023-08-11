source('./common_functions.R')

# List all the files within subdirectories
files <- list.files(path = "../Data/GSE89623_RAW", recursive = TRUE, full.names = TRUE, pattern = "*.txt.gz")

for (file in files){
    print(file)
    # Read the count matrix for all 
    counts = read.delim(file,row.names=1)
    # Reverse stranded
    counts = counts
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
colnames(counts_all) = gsub('_.*','',basename(files))
# Convert to TPM
#CPM to FPKM
counts_all = CPMToFPKM(counts_all)
#RPKM to TPM
counts_all = FPKMToTPM(counts_all)
# Save the counts matrix
exportCount(counts_all, 'GSE89623', 'GSE89623')
