library(Seurat)
library(trqwe)
library(dplyr)
library(data.table)
library(ggplot2)

# Read the GSE ID from command line
GSE = commandArgs(trailingOnly=TRUE)
# test if there is only one argument: if not, return an error
if (length(GSE)!=1) {
  stop("One argument must be supplied (GSE ID)", call.=FALSE)
} 

get_sig = function(sigs,suff){
    # If suff doesnt contain - then return the genes in the signature set
    if (!grepl('-', suff)){
        sigs = sigs[,grep(suff, colnames(sigs))]
    } else {
    # Split by - and return the genes in each signature
        sigs = lapply(strsplit(suff,'-')[[1]], function(x) sigs[,grep(x, colnames(sigs))])
    }
    return(sigs)
}

gbm_tsne = function(GSE, GSM, suff){
    print(paste(GSM, suff))
    # Read the counts
    counts = mcreadRDS(paste0("./Data_generated/", GSE, "/Counts/", GSM, "_counts.rds"), mc.cores=4)
    # Read the AUCell scores
    aucs = read.csv(paste0("./Output/", GSE, "/AUCell/", GSM, "-AUCell.csv"), row.names = 1)
    # Read the signatures
    sigs = read.csv("./Signatures/GBM_signatures.csv")
    #split suff by - and get the genes in each signature
    sigs = get_sig(sigs, suff)
    genes = unique(unlist(sigs, use.names = F))
    genes = genes[genes %in% rownames(counts)]
    # Subset only the signature genes
    counts = counts[genes,]
    # Select only the AUCell scores startin with suff
    aucs = aucs[,grep(suff, colnames(aucs))]
    # Set identity for each row by max AUCell score
    aucs = apply(aucs, 1, function(x) colnames(aucs)[which.max(x)])
    if (!all(dim(counts)) || is.null(dim(counts))){return()}
    # Convert to "dgCMatrix"
    counts = Matrix::Matrix(as.matrix(counts),sparse = T)
    # Make seurat object with normalized counts
    counts = CreateSeuratObject(counts, project = GSM)
    # Add Meta data of aucs
    counts = AddMetaData(counts, aucs, col.name = "Celltype")
    # Find variable genes
    counts = FindVariableFeatures(counts, selection.method = "vst", nfeatures = 2000, layer = "counts")
    # Scale the data
    counts <- ScaleData(counts)
    # Run PCA
    counts = RunPCA(counts, layer = "counts")
    # Run tSNE
    counts = RunTSNE(counts, dims = 1:2, check_duplicates = FALSE)
    # Run UMAP
    counts = RunUMAP(counts, dims = 1:2)
    # Define the colormapping for the AUCell scores
    # "NefNPC": #1f77b4,"NefAC" : #ff7f0e,"NefOPC" : #2ca02c, "NefMES" : #d62728
    cols=c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#1f77b4","#ff7f0e","#2ca02c","#d62728")
    names(cols) = c("NefNPC", "NefAC", "NefOPC", "NefMES","VerPN", "VerCL", "VerNL", "VerMES")
    # Plot the tSNE plot
    DimPlot(counts, reduction = "tsne", pt.size = 0.9, 
    group.by = "Celltype", cols = cols)
    ggsave(paste0("./figures/", GSE, "/PCA/", GSM, "_", suff, "_tsne.svg"), width = 10, height = 10)
    # Plot the UMAP plot
    DimPlot(counts, reduction = "umap", pt.size = 0.5, group.by = "Celltype", cols=cols)
    ggsave(paste0("./figures/", GSE, "/PCA/", GSM, "_", suff, "_umap.svg"), width = 10, height = 10)
}

# GSM file list
GSMs = list.files(paste0("./Data_generated/", GSE, "/Counts/"), "*_counts.rds")  %>% gsub('_counts.rds','',.)

# Create directory for output
dir.create(paste0("./Output/", GSE, "/PCA"), showWarnings = F, recursive = T)

# Iterate over GSM samples and generate rds
for (i in 1:length(GSMs)){
    gbm_tsne(GSE, GSMs[i], 'Nef')
    gbm_tsne(GSE, GSMs[i], 'Ver')
}
