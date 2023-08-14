source('./common_functions.R')
library(biomaRt)

# List all the files within subdirectories
files <- list.files(path = "../Data/GSE74557_RAW", recursive = TRUE, full.names = TRUE, pattern = "*.txt.gz")

for (file in files){
    print(file)
    # Read the count matrix for all 
    counts = read.delim(file, row.names=1)
    # Reverse stranded
    counts = counts[,5,drop=F]
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
colnames(counts_all) = gsub('_.*','',basename(files))
# Convert NCBI gene number to gene symbol
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genesymb = getBM(attributes = c("entrezgene_id","external_gene_name"),filters = 'entrezgene_id', values = rownames(df), mart = mart)
df = merge(df, genesymb, by.x = "row.names", by.y = "entrezgene_id")[,-1]
# merge duplicate gene names
df = df %>%
  group_by(external_gene_name) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as geneids
rownames(df) = df[,1]
df = df[,-1]
# Remove empty gene name
df = df[!rownames(df) == "",]
# Save the counts matrix
exportCount(counts_all, 'GSE74557', 'GSE74557')
