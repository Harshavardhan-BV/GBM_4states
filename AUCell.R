library(AUCell)
library(GSEABase)
library(trqwe)
library(data.table)
library(dplyr)

# Read the GSE ID from command line
GSE = commandArgs(trailingOnly=TRUE)
# test if there is only one argument: if not, return an error
if (length(GSE)!=1) {
  stop("One argument must be supplied (GSE ID)", call.=FALSE)
} 

run_AUCell <- function(GSE, GSM){
    print(GSM)
    # Read the counts
    matrix_AUCell = mcreadRDS(paste0("./Data_generated/",GSE,"/Counts/",GSM,"_counts.rds"), mc.cores=4) %>% as.matrix()
    # Get the gene sets
    geneSets = getGmt('./Signatures/GBM_states.gmt')
    # Subset only the expressed genes
    geneSets <- subsetGeneSets(geneSets, rownames(matrix_AUCell)) 
    # Build cell rankings
    matrix_AUCell <- AUCell_buildRankings(matrix_AUCell, plotStats = F)
    gc()
    # Calculate AUC metrics
    matrix_AUCell <- AUCell_calcAUC(geneSets, matrix_AUCell, nCores = parallel::detectCores()-2)
    # AUCmatrix has scores, row = pathway, col = cell
    matrix_AUCell = t(matrix_AUCell@assays@data$AUC) 
    matrix_AUCell = as.data.table(cbind(rownames(matrix_AUCell),matrix_AUCell))
    colnames(matrix_AUCell)[1] = 'Cellname'
    gc(reset = TRUE)
    if (any(dim(matrix_AUCell))){
        # Save the AUCell scores as a csv
        fwrite(matrix_AUCell, paste0("./Output/",GSE,"/AUCell/",GSM,"-AUCell.csv"))
    }
}

# GSM file list
GSMs = list.files(paste0('./Data_generated/',GSE,'/Counts/'),'*_counts.rds') %>% gsub('_counts.rds','',.)

# Create directory for output
dir.create(paste0("./Output/", GSE, "/AUCell"), showWarnings = F, recursive = T)

# Iterate over GSM samples and generate rds
for (GSM in GSMs){
    run_AUCell(GSE, GSM)
}
