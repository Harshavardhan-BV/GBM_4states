source('./common_functions.R')
library(biomaRt)

# List all the files within subdirectories
files <- list.files(path = "../Data/GSE100351_RAW", recursive = TRUE, full.names = TRUE, pattern = "*.RPKM.gz")

for (file in files){
    print(file)
    # Read the count matrix for all 
    counts = read.delim(file, row.names = 2)
    # preseve the rownames when subsetting RPKM
    counts = counts[,4,drop=FALSE]
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
# Select the Sample names from the file names
colnames(counts_all) = gsub("_.*", "", basename(files))
#TranscriptID to GeneName
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genesymb = getBM(attributes = c("ensembl_transcript_id_version","external_gene_name"),filters = 'ensembl_transcript_id_version', values = rownames(counts_all), mart = mart)
counts_all = merge(counts_all, genesymb, by.x = "row.names", by.y = "ensembl_transcript_id_version")[,-1]
# merge duplicate gene names
counts_all = counts_all %>%
  group_by(external_gene_name) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as geneids
rownames(counts_all) = counts_all[,1]
counts_all = counts_all[,-1]
# Convert to TPM
counts_all = FPKMToTPM(counts_all)
# Remove empty gene name
counts_all = counts_all[!rownames(counts_all)=='',]
# Save the counts matrix
exportCount(counts_all, 'GSE100351', 'GSE100351')
