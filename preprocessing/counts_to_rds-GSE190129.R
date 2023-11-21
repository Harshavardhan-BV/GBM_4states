source("./common_functions.R")
# List the files
files = list.files('../Data/GSE190129_RAW', full.names = T, pattern = '*.csv.gz')
for (file in files){
    # Get the GSM
    GSM = strsplit(basename(file), '_')[[1]][1]
    print(GSM)
    # Read the counts
    counts = read.csv(file, row.names=1)
    # Prepend GSM to the column names
    colnames(counts) = paste(GSM, colnames(counts), sep = '_')
    if (exists("counts_all")){
        counts_all = merge(counts_all, counts, by = "row.names", all = T)
        # Set rownames
        rownames(counts_all) = counts_all[,1]
        counts_all = counts_all[,-1]
        # Set nan to 0
        counts_all[is.na(counts_all)] = 0 
    } else {
        counts_all = counts
    }
}

# Perform QC and normalization
counts_all = SC_QC(counts_all)
# Save the counts matrix
exportCount(counts_all, 'GSE190129', 'GSE190129', sc=T)
