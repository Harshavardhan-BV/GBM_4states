source('./common_functions.R')

# List all the files within subdirectories
files <- list.files(path = "../Data/GSE195681_RAW", recursive = TRUE, full.names = TRUE, pattern = "*.out.tab.gz")

for (file in files){
    print(file)
    # Read the count matrix for all 
    counts = read.delim(file, header = F, skip = 4)
    # Reverse stranded
    counts = counts[,c(1,3)]
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
# Set the row names as geneids
rownames(counts_all) = counts_all[,1]
counts_all = counts_all[,-1]
# Select the Sample names from the file names
colnames(counts_all) = gsub('_ReadsPerGene.out.tab.gz','',basename(files))
# GeneID to GeneName
counts_all = GeneIDToSymb(counts_all)
# Convert to TPM
counts_all = countToTpm(counts_all)
# Save the counts matrix
exportCount(counts_all, 'GSE195681', 'GSE195681')
