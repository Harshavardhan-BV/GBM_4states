library(dplyr)
library(trqwe)
library(data.table)
# Read the count matrix for all 
counts = read.csv('../Data/CCLE/OmicsExpressionProteinCodingGenesTPMLogp1.csv')
rownames(counts) = counts[,1]
counts = counts[,-1]
# Read the names of cellines for GBM
celllines = read.csv('../Data/CCLE/cell lines in Glioblastoma.csv')
celllines = unlist(celllines['Depmap.Id'], use.names = F)
# Get the values of counts for index in celllines
counts <- counts[which(rownames(counts) %in% celllines),]
# removes the part after .. for columns
colnames(counts) <- gsub("\\.\\..*", "", colnames(counts))
# convert to dataframe
counts = data.frame(t(counts))
# Export the counts for the GBM celllines only as tsv
fwrite(counts, '../Data_generated/CCLE/Counts/CCLE_counts.tsv', sep='\t', row.names = T)
# Export the counts for GBM celllines as rds
mcsaveRDS(counts, "../Data_generated/CCLE/Counts/CCLE_counts.rds", mc.cores=4)
