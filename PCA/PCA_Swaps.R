library(trqwe)
library(dplyr)
library(data.table)

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

replic_pca = function(counts,sigs,houskeep,tgt,nswaps){
    # Select nswap genes from sigs
    nrow_sig = sigs[,tgt] %>% subset(., . != "") %>% length()
    rowid = sample(1:nrow_sig, nswaps, replace = F)
    # Select nswap genes from housekeeping
    hk = sample(houskeep, nswaps, replace = F)
    # Replace the genes in sigs with the genes in housekeeping
    sigs[rowid,tgt] = hk
    # Get the genes in the signature
    genes = unique(unlist(sigs, use.names = F))
    genes = genes[genes %in% rownames(counts)]
    # Subset only the signature genes
    counts = counts[genes,]
    if (!all(dim(counts)) || is.null(dim(counts))){return()}
    counts = counts[!apply(counts,1,sd)==0,]
    # Run PCA
    pca = prcomp(t(counts), center = T, scale. = T)
    # Return the explained variance of PC1
    exp_var = (pca$sdev^2 / sum(pca$sdev^2))[1:4]
    return (exp_var)
}

swap_pca = function(GSE, GSM, suff, nrepl=10){
    print(paste(GSM, suff))
    # Read the counts
    counts = mcreadRDS(paste0("../Data_generated/", GSE, "/Counts/", GSM, "_counts.rds"), mc.cores=4)
    # Read the housekeeping genes
    houskeep = read.delim('../Signatures/Housekeeping_genes.txt')[,1] %>% unlist(n)
    # Read the signatures
    sigs = read.csv("../Signatures/GBM_signatures.csv")
    # split suff by - and get the genes in each signature
    sigs = get_sig(sigs, suff)
    mega_df = NULL
    # Iterate over the signatures and nswaps
    for (tgt in colnames(sigs)){
        for (pct in seq(0,100,10)){
            # Get the number of genes to swap
            nrow_sig = sigs[,tgt] %>% subset(., . != "") %>% length()
            nswaps = round(pct/100 * nrow_sig)
            # Run the replicates
            exp_var = replicate(nrepl, replic_pca(counts,sigs,houskeep,tgt, nswaps))
            # Appends to list
            mega_df = cbind(tgt, nswaps, pct, t(exp_var)) %>% rbind(mega_df)
        }
    }
    # Conver to data frame
    mega_df = as.data.frame(mega_df)
    colnames(mega_df) = c('Signature', 'Swaps','Swap_pct' , 'PC1_Var', 'PC2_Var', 'PC3_Var', 'PC4_Var')
    # Save the results
    fwrite(mega_df, paste0("../Output/", GSE, "/PCA_Swaps/", GSM, "_", suff, "_swaps.csv")) 
}

# GSM file list
GSMs = list.files(paste0("../Data_generated/", GSE, "/Counts/"), "*_counts.rds")  %>% gsub('_counts.rds','',.)

# Create directory for output
dir.create(paste0("../Output/", GSE, "/PCA_Swaps"), showWarnings = F, recursive = T)

# Iterate over GSM samples and generate rds
for (GSM in GSMs){
    swap_pca(GSE, GSM, 'Nef')
    swap_pca(GSE, GSM, 'Ver')
}