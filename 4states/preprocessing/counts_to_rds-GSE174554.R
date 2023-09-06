source("./common_functions.R")
# List the files
mtxs = list.files('../Data/GSE174554_RAW', full.names = T, pattern = '*.mtx.gz')
barcodes =  list.files('../Data/GSE174554_RAW', full.names = T, pattern = '*barcodes.tsv.gz')
features = list.files('../Data/GSE174554_RAW', full.names = T, pattern = '*features.tsv.gz')
# Get order for sorting
mtx_ord = order(gsub('_matrix.mtx.gz','', basename(mtxs)))
bcd_ord = order(gsub('_barcodes.tsv.gz','', basename(barcodes)))
ftr_ord = order(gsub('_features.tsv.gz','', basename(features)))
# Sort the values
mtxs = mtxs[mtx_ord]
barcodes = barcodes[bcd_ord]
features = features[ftr_ord]

for (i in 1:length(mtxs)){
    # Get the GSM
    GSM = strsplit(basename(mtxs[i]), '_')[[1]][1]
    print(GSM)
    # Read the counts
    counts = ReadMtx(mtx = mtxs[i], cells = barcodes[i], features = features[i], feature.column = 1)
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
# Perform QC and normalization
counts_all = SC_QC(counts_all)
# Save the counts matrix
exportCount(counts_all, 'GSE174554', 'GSE174554', sc=T)
