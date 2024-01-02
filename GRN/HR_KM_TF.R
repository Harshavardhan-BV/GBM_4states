library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(trqwe)
library(grid)
library(ggthemes)

theme_Publication <- function(base_size=14, base_family="sans") {
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.2), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = "#ffffff"),
                plot.background = element_rect(colour = "#ffffff"),
                panel.border = element_rect(colour = "black"),
                axis.title = element_text(face = "bold",size = rel(1.5)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(size = rel(1.2)),
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
}

# Read the GSE ID from command line
GSE = commandArgs(trailingOnly=TRUE)
# test if there is only one argument: if not, return an error
if (length(GSE)!=1) {
  stop("One argument must be supplied (GSE ID)", call.=FALSE)
} 

# add method to grid.draw
grid.draw.ggsurvplot <- function(x){
survminer:::print.ggsurvplot(x, newpage = FALSE)
}

split_high_low = function(x) {
  ifelse(x > median(x, na.rm = TRUE), "High", "Low") %>% as.factor() %>% relevel(x, ref='Low')
}

# Plot the survival curves
KM_plot = function(survdat, gene, survtype, GSE, sig){
    # Bullshit R. GARBAGE ASS HORRIBLE PIECE OF SOFTWARE 
    survtype <<- survtype
    gene <<- gene
    # Drop the NAs of survtype
    survdat = survdat[!is.na(survdat[,survtype]),]
    # Split to high low based on median
    survdat[,gene] = split_high_low(survdat[,gene])
    # Fit the model to data limit
    print(paste(survtype, gene))
    fit <- survfit(as.formula(paste0("Surv(", survtype,".time,", survtype,") ~ ", gene)), data = survdat)
    res.cox <- coxph(as.formula(paste0("Surv(", survtype,".time,", survtype,") ~ ", gene)), data = survdat)
    res.cox = summary(res.cox)$coefficients[,2]
    max_time = max(survdat[,paste0(survtype,".time")]) * 0.8
    # Make a KM plot
    splots = ggsurvplot(fit, data = survdat, pval = T, ggtheme = theme_Publication(base_size = 16), conf.int = T, pval.coord = c(max_time,0.95), pval.size=4)
    # Add the HR to plot
    splots$plot = splots$plot + annotate("text", x=max_time, y=1, label= paste0('Hazard Ratio = ',format(res.cox, digits=3)))
    # Save the figure
    ggsave(paste0("../figures/", GSE, "/GRN/survival/",gene,'_',survtype,'_KMplot.svg'), splots, width = 8, height = 7)
}

HR_plot = function(survdat, genes, survtype, GSE, GSM, sig){
    # Drop the NAs of survtype
    survdat = survdat[!is.na(survdat[,survtype]),]
    # Split into high or low based on median for each gene
    survdat[, genes$V1] <- lapply(survdat[, genes$V1], split_high_low)
    # Get the hazard ratio for each gene
    df_res = list()
    for (i in 1:length(genes$V1)){
        # Do the cox regression
        gene = genes[i,1] 
        res.cox <- coxph(as.formula(paste0("Surv(", survtype,".time,", survtype,") ~ ", gene)), data = survdat)
        df_res = rbind(df_res, c(gene, summary(res.cox)$coefficients, genes[i,2]))
    }
    df_res = data.frame(df_res)
    colnames(df_res) = c('Gene', 'Coef', 'Exp_Coef', 'SE', 'z', 'pvalue', 'Subtype')
    # Convert to tibble and make the columns 2:-2 numeric
    df_res = df_res %>% as_tibble() %>% mutate_at(3:ncol(df_res)-1, as.numeric) %>% mutate_at(c(1, ncol(df_res)), as.character)
    # Make a HR plot
    HRplot = ggforestplot::forestplot(df_res, name=Gene, estimate = Coef, se=SE, pvalue= pvalue, psignif = 0.05, colour = Subtype, logodds = T, xlab = 'Hazard Ratio') + scale_color_manual(values = c("NefMES" = "#d62728","NefNPC"="#1f77b4","VerMES"="#d62728", "VerPN"="#1f77b4")) + theme_Publication()
    # Bleh
    ggsave(paste0('../figures/', GSE,'/GRN/survival/',GSM,'_',sig,'_',survtype,'_HRplot.svg'), width = 7, height = 7.5)
}

survival_analysis = function(GSE, GSM, sig){
    # Read the survival data
    survdat = read.delim('../Data/TCGA/GBM_survival.txt', sep='\t', row.names=1)
    # Read the TCGA expression data
    counts = mcreadRDS('../Data_generated/TCGA/Counts/TCGA_counts.rds', mc.cores=4) %>% t() %>% data.frame()
    # Read the top gene signatures
    files =  list.files(paste0('../Output/',GSE,'/GRN/'),pattern = paste0(GSM,'_',sig,'.*-topTFs.tsv'), full.names = TRUE)
    # Make a dataframe of TFs & the subtype
    for (flies in files){
        TF1 = read.delim(flies, sep='\t', header = F)
        TF1$subtype = basename(flies) %>% gsub('-topTFs.tsv','',.) %>% gsub(paste0(GSM,'_'),'',.)
        if (exists('TFs')){
            TFs = rbind(TFs, TF1)
        } else {
            TFs = TF1
        }
    }
    # Select the genes from the counts
    counts = counts[,TFs$V1]
    # Replace the . with - in the rownames
    rownames(counts) = gsub('\\.', '-', rownames(counts))
    # Select the string that matches the regex pattern
    counts$ID = gsub('(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[A-Z0-9]{2}).*', '\\1', rownames(counts))
    # Merge the ssGSEA scores with the survival data
    survdat = merge(survdat, counts, by.x = "row.names", by.y = 'ID', all.y = T)
    # Remove rows with NA in Row.names
    survdat = survdat[!is.na(survdat$Row.names),]
    # Do HR plot of all the genes
    HR_plot(survdat, TFs, 'OS', GSE, GSM, sig)
    HR_plot(survdat, TFs, 'PFI', GSE, GSM, sig)
    # Do KM plot for each gene
    for (gene in TFs$V1){
        KM_plot(survdat, gene, 'OS', GSE)
    }
}


GSMs = list.files(paste0('../Data_generated/',GSE,'/Counts/'),'*_counts.rds') %>% gsub('_counts.rds','',.)

# Create directory for output
dir.create(paste0("../figures/", GSE, "/GRN/survival"), showWarnings = F, recursive = T)

# Iterate over GSM samples
for (GSM in GSMs){
    survival_analysis(GSE, GSM, 'Nef')
    survival_analysis(GSE, GSM, 'Ver')
}
