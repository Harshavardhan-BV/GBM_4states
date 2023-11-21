source("./common_functions.R")
# List the files
files = list.files('../Data/GSE174470_RAW', full.names = T, pattern = '*.tsv.gz')

for (file in files){
    # Get the GSM
    GSM = strsplit(basename(file), '_')[[1]][2] %>% gsub('.tsv.gz','',.)
    print(GSM)
    # Read the counts
    counts = read.delim(file, row.names = 1)
    # Prepend GSM to the column names
    colnames(counts) = paste(GSM, colnames(counts), sep = '_')
    # Convert gene IDs to symbols
    counts = GeneIDToSymb(counts)
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
# Remove genes with zero counts
counts_all = counts_all[,colSums(counts_all)>0]
# Perform QC and normalization
counts_all = SC_QC(counts_all)
# to upper
rownames(counts_all) = toupper(rownames(counts_all))
# Save the counts matrix
exportCount(counts_all, 'GSE174470', 'GSE174470', sc=T)
