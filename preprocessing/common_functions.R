library(dplyr)
library(trqwe)
library(data.table)

# Save the processed counts matrix
exportCount <- function(counts, GSE, GSM, sc=FALSE){
    # Create the directory for the counts
    dir.create(paste0("../Data_generated/", GSE, "/Counts"), showWarnings = FALSE, recursive = TRUE)
    # Save the counts matrix as rds
    mcsaveRDS(counts, paste0("../Data_generated/", GSE, "/Counts/", GSM, "_counts.rds"), mc.cores=4)
    if (!sc){
        # If bulk save as tsv as well
        fwrite(counts, paste0("../Data_generated/", GSE, "/Counts/", GSM, "_counts.tsv"), sep="\t", row.names = T)
    }
}
# GeneID to GeneSymbol
GeneIDToSymb <- function(counts){
    # Remove the version number(?) from the gene IDs
    geneid = gsub("\\..*", "", rownames(counts))
    # Read the gene annotation file
    hg_gene_names = read.table(file="../Signatures/hg38_gene_length_kb.txt", sep = "\t",header = F,stringsAsFactors = F)
    mm_gene_names = read.table(file="../Signatures/mm10_gene_length_kb.txt", sep = "\t",header = F,stringsAsFactors = F)
    # Check if the counts are from human or mouse
    if (length(intersect(geneid, hg_gene_names$V1)) > 0){
        gene_names = hg_gene_names
    } else if (length(intersect(geneid, mm_gene_names$V1)) > 0){
        gene_names = mm_gene_names
    } else {
        print("Error: Gene names not found")
    }
    # Create a new column with gene names according to the mapping on the gene length file
    counts$GeneSymbol = gene_names[match(geneid, gene_names$V1), 2]
    # Remove nan gene names
    counts = counts[!is.na(counts$GeneSymbol),]
    # Merge duplicate row names and take mean of numeric values
    counts = counts %>%
    group_by(GeneSymbol) %>%
    summarise_all(mean) %>% data.frame()
    # Set the row names as gene names
    rownames(counts) = counts$GeneSymbol
    # Remove the GeneSymbol column
    counts = counts[,-1]
    # Return the counts matrix
    return(counts)
}

# TPM
countToTpm <- function(rawCount) {
    rawCount <- cbind(rownames(rawCount), rawCount)
    # Read gene length annotation
    geneLengthMap_hg = read.table(file="../Signatures/hg38_gene_length_kb.txt", sep = "\t",header = F,stringsAsFactors = F)
    commonGenes_hg = intersect(rawCount[,1], geneLengthMap_hg[,2])
    geneLengthMap_mm = read.table(file="../Signatures/mm10_gene_length_kb.txt", sep = "\t",header = F,stringsAsFactors = F)
    commonGenes_mm = intersect(rawCount[,1], geneLengthMap_mm[,2])
    geneCol = 2
    # Check which species has more common genes and use that for TPM calculation
    if(length(commonGenes_hg) > length(commonGenes_mm)) {
        geneLengthMap = geneLengthMap_hg
        commonGenes = commonGenes_hg
    }
    else{
        geneLengthMap = geneLengthMap_mm
        commonGenes = commonGenes_mm
    }
    # Get the common gene indices for count matrix
    idx1 = match(commonGenes,rawCount[,1])
    # Subset count matrix for common genes
    rawCount = rawCount[idx1,-c(1)]
    # Get the common gene indices for gene length annotation
    idx2 = match(commonGenes,geneLengthMap[,2])
    # Subset gene length annotation for common genes
    geneLengthSubset = geneLengthMap[idx2, ]
    geneSymbol = geneLengthSubset[ ,geneCol]
    featureLength = geneLengthSubset[ ,3]

    rownames(rawCount) = geneSymbol
    CellNames = colnames(rawCount) 
    # Ensure valid arguments.
    stopifnot(length(featureLength) == nrow(rawCount))

    # Compute effective lengths of features in each library.
    effLen <- featureLength

    # Process one column at a time.
    tpm <- do.call(cbind, lapply(1:ncol(rawCount), function(i) {
    rate = (rawCount[,i])/effLen
    rate/sum(rate, na.rm = T) * 1e6
    }))

    rm(rawCount)
    gc(verbose = F)
    # Convert to log2tpm
    tpm = log2(tpm+1)
    gc(verbose = F)
    # Copy the row and column names from the original matrix.
    colnames(tpm) <- CellNames 
    rownames(tpm) <- toupper(geneSymbol)
    # Set na values to 0
    tpm[is.na(tpm)] = 0
    # Return the TPM matrix
    return(data.frame(tpm))
}

FPKMToTPM <- function(counts) {
    # Store the rownames and colnames
    geneSymbol = rownames(counts)
    CellNames = colnames(counts)
    # Process one column at a time.
    tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = sum(counts[,i])
    (counts[,i]/rate) * 1e6
    }))
    rm(counts)
    gc(verbose = F)
    # Convert to log2tpm
    tpm = log2(tpm+1)
    gc(verbose = F)
    # Set the row and column names
    colnames(tpm) <- CellNames 
    rownames(tpm) <- geneSymbol
    # Return the TPM matrix
    return(data.frame(tpm))
}