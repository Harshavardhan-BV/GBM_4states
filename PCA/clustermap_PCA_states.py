import os
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# Read the GSE ID from command line
parser = argparse.ArgumentParser(
                    prog='clustermap_PCA_states.py',
                    description='Plots a clustermap of the PC1,PC2 loadings of the PCA')
parser.add_argument('GSE', type=str, nargs=1, help='GSE ID')
GSE = parser.parse_args().GSE[0]

def loada_clustermap(GSE, GSM, suff):
    print(GSM,suff)
    # Load the loadings matrix
    load_df = pd.read_csv('../Output/'+GSE+'/PCA/'+GSM+'_loadings_'+suff+'.tsv',sep='\t',index_col=0)
    # Load the variance explained
    var_df = pd.read_csv('../Output/'+GSE+'/PCA/'+GSM+'_expvar_'+suff+'.tsv',sep='\t',header=None)
    var_df = var_df*100
    gbmgenes = pd.read_csv("../Signatures/GBM_signatures.csv")
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
    plt.savefig('../figures/'+GSE+'/PCA/'+GSM+'_loadings_clustermap_'+suff+'.svg')
    plt.close()
    plt.clf()

# Make the output directory
os.makedirs('../figures/'+GSE+'/PCA/', exist_ok=True)
# Iterative over the signatures
for suff in ['Nef','Ver']:
    # List files in Output
    files = glob.glob('../Output/'+GSE+'/PCA/*_loadings_'+suff+'.tsv')
    # Get the GSMs
    GSMs = [os.path.basename(x).replace(f'_loadings_'+suff+'.tsv','') for x in files]
    # Iterate over GSM samples
    for GSM in GSMs:
        loada_clustermap(GSE,GSM,suff)