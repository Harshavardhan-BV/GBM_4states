source("./common_functions.R")
# List the files
files = list.files('../Data/GSE200984_RAW', full.names = T, pattern = '*.txt.gz')
for (file in files){
    # Get the GSM
    GSM = strsplit(basename(file), '_')[[1]][1]
    print(GSM)
    # Read the counts
    tryCatch({
        # Read the counts
        counts = read.table(file, row.names = 1)
    }, error = function(e) {
        counts = read.table(file, row.names = NULL)
        # Merge duplicate gene names
        counts = counts %>%
        group_by(counts[,1]) %>%
        summarize_all(mean) %>% data.frame()
        # And set it as rownames
        rownames(counts) = counts[,1]
        counts = counts[,-1]
    }
    )

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
counts = SC_QC(counts)
# Save the counts matrix
exportCount(counts, 'GSE200984', 'GSE200984', sc=T)
