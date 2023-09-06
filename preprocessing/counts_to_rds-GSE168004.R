source('./common_functions.R')

# Cell lines WT
files = c('mgg18','mgg23','mgg75','mgh143')
for (file in files){
    # Read the counts
    counts = read.csv(paste0("../Data/GSE168004_RAW/GSE168004_",file,"_tpm.csv.gz"),row.names = 1)
    if (exists("counts_all")){
        counts_all = merge(counts_all, counts, by = 'row.names', all = T)
        rownames(counts_all) = counts_all[,1]
        counts_all = counts_all[,-1]
        counts_all[is.na(counts_all)] = 0
    } else {
        counts_all = counts
    }
}
# Perform QC
counts_all = SC_QC(counts_all, normalize = F)
# Save the counts matrix
exportCount(counts_all, 'GSE168004', 'celllines', sc=T)

# Mouse
counts = read.csv("../Data/GSE168004_RAW/GSE168004_mouse_tpm.csv.gz",row.names = 1)
# Perform QC
counts = SC_QC(counts, normalize = F)
# Convert to uppercase
rownames(counts) = toupper(rownames(counts))
# Save the counts matrix
exportCount(counts, 'GSE168004', 'mouse', sc=T)

# Cell lines OSM
counts = read.csv("../Data/GSE168004_RAW/GSE168004_OSM_celllines_tpm.csv.gz",row.names = 1)
# Perform QC
counts = SC_QC(counts, normalize = F)
# Save the counts matrix
exportCount(counts, 'GSE168004', 'OSM_celllines', sc=T)

# Mouse OSM
counts = read.csv("../Data/GSE168004_RAW/GSE168004_OSM_mouse_tpm.csv.gz",row.names = 1)
# Perform QC
counts = SC_QC(counts, normalize = F)
# Convert to uppercase
rownames(counts) = toupper(rownames(counts))
# Save the counts matrix
exportCount(counts, 'GSE168004', 'OSM_mouse', sc=T)
