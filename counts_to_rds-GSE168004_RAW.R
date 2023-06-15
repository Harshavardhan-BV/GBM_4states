library(dplyr)
library(trqwe)

df = read.csv("./Data/GSE168004_RAW/GSE168004_OSM_celllines_tpm.csv.gz",row.names = 1)
mcsaveRDS(df, "./Data_generated/GSE168004/OSM_celllines_counts.rds", mc.cores=4)

df = read.csv("./Data/GSE168004_RAW/GSE168004_mgg23_tpm.csv.gz",row.names = 1)
mcsaveRDS(df, "./Data_generated/GSE168004/mgg23_counts.rds", mc.cores=4)

df = read.csv("./Data/GSE168004_RAW/GSE168004_mgg75_tpm.csv.gz",row.names = 1)
mcsaveRDS(df, "./Data_generated/GSE168004/mgg75_counts.rds", mc.cores=4)

