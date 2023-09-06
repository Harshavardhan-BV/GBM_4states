source("./common_functions.R")
# List the files
mtxs = list.files('../Data/GSE182109_RAW', full.names = T, pattern = '*.mtx.gz')
barcodes =  list.files('../Data/GSE182109_RAW', full.names = T, pattern = '*barcodes.tsv.gz')
features = list.files('../Data/GSE182109_RAW', full.names = T, pattern = '*features.tsv.gz')
# Get order for sorting
mtx_ord = order(gsub('_matrix.mtx.gz','', basename(mtxs)))
bcd_ord = order(gsub('_barcodes.tsv.gz','', basename(barcodes)))
ftr_ord = order(gsub('_features.tsv.gz','', basename(features)))
# Sort the values
mtxs = mtxs[mtx_ord]
barcodes = barcodes[bcd_ord]
features = features[ftr_ord]
# Read metadata
metdat = read.delim('../Data/GSE182109_RAW/Meta_Data_GBMatlas.txt')
# Select only glioma cells
glioma = metdat[metdat$Assignment=='Glioma',]
for (i in 1:length(mtxs)){
    # Get the GSM
    GSM = strsplit(basename(mtxs[i]), '_')[[1]][1]
    print(GSM)
    # Read the counts
    counts = ReadMtx(mtx = mtxs[i], cells = barcodes[i], features = features[i])
    # Prepend GSM to the column names
    colnames(counts) = paste(GSM, colnames(counts), sep = '_')
    # Select only glioma cells
    counts = counts[,colnames(counts) %in% rownames(glioma)]
    if (exists("counts_all")){
        counts_all = merge(counts_all, counts, by = 'row.names', all = T)
        rownames(counts_all) = counts_all[,1]
        counts_all = counts_all[,-1]
        counts_all[is.na(counts_all)] = 0
    } else {
        counts_all = counts
    }
}
# Perform QC and normalization
counts_all = SC_QC(counts_all)
# Save the counts matrix
exportCount(counts_all, 'GSE182109', 'GSE182109', sc=T)
