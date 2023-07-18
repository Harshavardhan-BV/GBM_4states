import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''

def loada_clustermap(GSE, GSM, suff):
    print(GSM,suff)
    # Load the loadings matrix
    load_df = pd.read_csv('./Output/'+GSE+'/PCA/'+GSM+'_loadings_'+suff+'.tsv',sep='\t',index_col=0)
    # Load the variance explained
    var_df = pd.read_csv('./Output/'+GSE+'/PCA/'+GSM+'_expvar_'+suff+'.tsv',sep='\t',header=None)
    var_df = var_df*100
    gbmgenes = pd.read_csv("./Signatures/GBM_signatures.csv")
    # Select only the columns that start with suff
    gbmgenes = gbmgenes.filter(like=suff)
    # select only the genes in the signature
    genes = gbmgenes.melt().dropna()
    genes = genes[genes.value.isin(load_df.index)]
    # make a colour map for the genes
    lut = { col: f'tab:{clr}' for col, clr in zip(gbmgenes.columns,['red','orange','green','blue'])}
    # map the colours to the genes
    row_colors = load_df.index.map(genes.drop_duplicates(subset=['value']).set_index('value')['variable'])
    row_colors = row_colors.map(lut)
    row_colors = row_colors.fillna('white')
    # Plot the loadings matrix
    g = sns.clustermap(data=load_df, cmap='coolwarm',center=0, row_colors=row_colors, col_cluster=False)
    # create a list of patches for the legend
    patches = [
        mpatches.Patch(color=color, label=column) for column, color in lut.items()
    ]
    # add legend for row_colours given by lut
    plt.legend(handles = patches , bbox_to_anchor=(2, 1), loc='upper left')
    # Append the variances explained to x-ticks
    xticklabels = [x+'\n'+str(round(y[0],2))+'%' for x,y in zip(load_df.columns,var_df.values)]
    # Set the x-ticks
    g.ax_heatmap.set_xticklabels(xticklabels,rotation=0)
    plt.savefig('./figures/'+GSE+'/PCA/'+GSM+'_loadings_clustermap_'+suff+'.svg')
    plt.close()
    plt.clf()

# Sort the loadings of PC1 and plot a barplot. Colour the bars according to the geneset
def loada_barplot(GSE, GSM, suff):
    print(GSM,suff)
    # Load the loadings matrix
    load_df = pd.read_csv('./Output/'+GSE+'/PCA/'+GSM+'_loadings_'+suff+'.tsv',sep='\t',index_col=0)
    # Load the variance explained
    var_df = pd.read_csv('./Output/'+GSE+'/PCA/'+GSM+'_expvar_'+suff+'.tsv',sep='\t',header=None)
    var_df = var_df*100
    # Read the GBM signature
    gbmgenes = pd.read_csv("./Signatures/GBM_signatures.csv")
    # Select only the columns that start with suff
    gbmgenes = gbmgenes.filter(like=suff)
    # select only the genes in the signature and merge with load_df
    gbmgenes = gbmgenes.melt(var_name='Signature').dropna()
    load_df = load_df.merge(gbmgenes,left_index=True,right_on='value', how='left')
    # Plot the loadings matrix but dont plot the ylabels
    g = sns.barplot(data=load_df,x='PC1',y=load_df.index, hue='Signature', dodge=False, orient='h', order= load_df.sort_values('PC1',ascending=False).index)
    g.set(yticklabels=[], yticks=[])
    plt.xlabel('PC1\n'+str(round(var_df.values[0][0],2))+'%')
    plt.tight_layout()
    plt.savefig('./figures/'+GSE+'/PCA/'+GSM+'_loadings_barplot_'+suff+'.svg')
    plt.close()
    plt.clf()

# GSE ID
GSE = "GSE168004"
# GSE = "GSE131928"
# GSE = "GSE182109"
# GSE = "CCLE"
# GSE = "TCGA"
# Make the output directory
os.makedirs('figures/'+GSE+'/PCA/', exist_ok=True)
# Iterative over the signatures
for suff in ['Nef','Ver']:
    # List files in Output
    files = glob.glob('Output/'+GSE+'/PCA/*_loadings_'+suff+'.tsv')
    # Get the GSMs
    GSMs = [os.path.basename(x).replace(f'_loadings_'+suff+'.tsv','') for x in files]
    # Iterate over GSM samples
    for GSM in GSMs:
        loada_clustermap(GSE,GSM, suff)
        loada_barplot(GSE,GSM, suff)
