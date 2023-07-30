source('./common_functions.R')
# List files in the directory
files = list.files(path = "../Data/GSE174470_RAW/", pattern = "*.tsv.gz")

for (file in files){
    print(file)
    # Read the counts
    df = read.delim(paste0("../Data/GSE174470_RAW/",file),row.names = 1)
    # GeneID to GeneName
    df = GeneIDToSymb(df)
    # Convert to TPM
    df = countToTpm(df)
    # uppercase gene names
    rownames(df) = toupper(rownames(df))
    # Save the counts matrix
    GSM = file %>% gsub('.tsv.gz','',.) %>% gsub('GSE174470_','',.)
    exportCount(df, 'GSE174470', GSM, sc=T)
}
