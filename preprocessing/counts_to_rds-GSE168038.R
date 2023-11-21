source("./common_functions.R")
# List the files
mtxs = list.files('../Data/GSE168038_RAW', full.names = T, pattern = '*.mtx.gz')
barcodes =  list.files('../Data/GSE168038_RAW', full.names = T, pattern = '*barcodes.tsv.gz')
features = list.files('../Data/GSE168038_RAW', full.names = T, pattern = '*features.tsv.gz')
for (i in 1:length(mtxs)){
    # Get the GSM
    GSM = strsplit(basename(mtxs[i]), '_')[[1]][1]
    print(GSM)
    # Read the counts
    counts = ReadMtx(mtx = mtxs[i], cells = barcodes[i], features = features[i])
    # Prepend GSM to the column names
    colnames(counts) = paste(GSM, colnames(counts), sep = '_')
    if (exists("counts_all")){
        if (nrow(counts_all)==nrow(counts)){
            counts_all = cbind(counts_all, counts)
        } else {
            print('Row number mismatch')
        }
    } else {
        counts_all = counts
    }
}
# Uppercase 
rownames(counts_all) = toupper(rownames(counts_all))
# Perform QC and normalization
counts_all = SC_QC(counts_all)
# Save the counts matrix
exportCount(counts_all, 'GSE168038', 'GSE168038', sc=T)
