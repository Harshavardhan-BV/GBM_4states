source("./common_functions.R")
# List the files
files = list.files('../Data/GSE130119_RAW', full.names = T, pattern = '*.csv.gz')
for (file in files){
    # Get the GSM
    GSM = strsplit(basename(file), '_')[[1]][1]
    print(GSM)
    # Read the counts
    counts = read.csv(file, row.names=1)
    # Prepend GSM to the column names
    colnames(counts) = paste(GSM, colnames(counts), sep = '_')
    # Convert GeneID to GeneSymb
    counts = GeneIDToSymb(counts) 
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

# Perform QC and normalization
counts_all = SC_QC(counts_all)
# Save the counts matrix
exportCount(counts_all, 'GSE130119', 'GSE130119', sc=T)
