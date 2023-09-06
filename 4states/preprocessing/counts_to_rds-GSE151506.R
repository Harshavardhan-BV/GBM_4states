source("./common_functions.R")
# List the files
files = list.files('../Data/GSE151506_RAW', full.names = T, pattern = 'GSM')
for (file in files){
    # Get the GSM
    GSM = strsplit(basename(file), '_')[[1]][1]
    print(GSM)
    # Read the counts
    counts = read.delim(file,row.names=1)
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
counts_all = SC_QC(counts_all,normalize = F)
# Save the counts matrix
exportCount(counts_all, 'GSE151506', 'GSM', sc=T)

rm(counts_all)

# Seperate files
files = list.files('../Data/GSE151506_RAW', full.names = T, pattern = 'GSE')
for (file in files){
    # Get the GSM
    GSM = strsplit(basename(file), '_')[[1]][2]
    print(GSM)
    # Read the counts
    counts = read.delim(file,row.names=1)
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
counts_all = SC_QC(counts_all,normalize = F)
# Save the counts matrix
exportCount(counts_all, 'GSE151506', 'GSE', sc=T)