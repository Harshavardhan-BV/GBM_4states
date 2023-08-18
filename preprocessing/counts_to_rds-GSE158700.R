source('./common_functions.R')
library(biomaRt)

# List all the files within subdirectories
files <- list.files(path = "../Data/GSE158700_RAW", recursive = TRUE, full.names = TRUE, pattern = "*.sf.txt.gz")

for (file in files[1:28]){
    print(file)
    # Read the count matrix for all 
    counts = read.delim(file, row.names=1)
    # Only counts
    counts = counts[,3,drop=F]
    # merge with the counts matrix
    if (exists("counts_all")){
        if (nrow(counts)==nrow(counts_all)){
            counts_all = cbind(counts_all,counts)
        } else {
            print(paste("Error: Row number mismatch", nrow(counts), nrow(counts_all)))
        }
    } else {
        counts_all = counts
    }
}
# Set the column names as sample names
colnames(counts_all) = gsub('_.*','',basename(files[1:28]))
# Transcript ID version fucking things up
rownames(counts_all) = gsub("\\..*", "", rownames(counts_all))
# Change geneid to genename
mart_mm = useMart("ensembl", dataset="mmusculus_gene_ensembl")
genesymb = getBM(attributes = c("ensembl_transcript_id","external_gene_name"),filters = 'ensembl_transcript_id', values = rownames(counts_all), mart = mart_mm)
counts_all = merge(counts_all, genesymb, by.x = "row.names", by.y = "ensembl_transcript_id")[,-1]
# merge duplicate gene names
counts_all = counts_all %>%
  group_by(external_gene_name) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as geneids
rownames(counts_all) = counts_all[,1]
counts_all = counts_all[,-1]
# Remove empty gene name
counts_all = counts_all[!rownames(counts_all) == "",]
# uppercase
rownames(counts_all) = toupper(rownames(counts_all))
# Save the counts matrix
exportCount(counts_all, 'GSE158700-B', 'mouse')

# There are human samples seperately
rm(counts_all)
for (file in files[29:length(files)]){
    print(file)
    # Read the count matrix for all 
    counts = read.delim(file, row.names=1)
    # Only counts
    counts = counts[,3,drop=F]
    # merge with the counts matrix
    if (exists("counts_all")){
        if (nrow(counts)==nrow(counts_all)){
            counts_all = cbind(counts_all,counts)
        } else {
            print(paste("Error: Row number mismatch", nrow(counts), nrow(counts_all)))
        }
    } else {
        counts_all = counts
    }
}
# Set the column names as sample names
colnames(counts_all) = gsub('_.*','',basename(files[29:length(files)]))
# Change geneid to genename
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
genesymb = getBM(attributes = c("ensembl_transcript_id_version","external_gene_name"),filters = 'ensembl_transcript_id_version', values = rownames(counts_all), mart = mart)
counts_all = merge(counts_all, genesymb, by.x = "row.names", by.y = "ensembl_transcript_id_version")[,-1]
# merge duplicate gene names
counts_all = counts_all %>%
  group_by(external_gene_name) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as geneids
rownames(counts_all) = counts_all[,1]
counts_all = counts_all[,-1]
# Remove empty gene name
counts_all = counts_all[!rownames(counts_all) == "",]
# Save the counts matrix
exportCount(counts_all, 'GSE158700-B', 'human')