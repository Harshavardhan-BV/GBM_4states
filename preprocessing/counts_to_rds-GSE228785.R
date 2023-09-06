source("./common_functions.R")
# List the files
files = list.files('../Data/GSE228785_RAW', full.names = T, pattern = '*.h5')
for (file in files){
    # Get the GSM
    GSM = strsplit(basename(file), '_')[[1]][1]
    print(GSM)
    # Read the counts
    counts = Read10X_h5(file)
    # Prepend GSM to the column names
    colnames(counts) = paste(GSM, colnames(counts), sep = '_')
    if (exists("counts_all")){
        if(nrow(counts_all)==nrow(counts)){
            counts_all = cbind(counts_all,counts)
        } else {
            print('Row number mismatch')
        }
    } else {
        counts_all = counts
    }
}
# Convert to uppercase
rownames(counts_all) = toupper(rownames(counts_all))
# Perform QC and normalization
counts_all = SC_QC(counts_all)
# Save the counts matrix
exportCount(counts_all, 'GSE228785', 'GSE228785', sc=T)
