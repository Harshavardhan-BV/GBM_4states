import os
import glob
import argparse
import itertools as it
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''

# Read the GSE ID from command line
parser = argparse.ArgumentParser(
                    prog='hmap_Expr.py',
                    description='Plots heatmaps of the correlation between expression of genes in signatures')
parser.add_argument('GSE', type=str, nargs=1, help='GSE ID')
GSE = parser.parse_args().GSE[0]

def corr_pair_plot(corr_df1, gbmgenes1, GSE, GSM):
    combinats = list(it.combinations(gbmgenes1.columns,2))
    for combi in combinats:
        suff = '-'.join(combi) 
        # Select the genes of geneset combination
        gbmgenes = gbmgenes1.loc[:,combi]
        # Select only the genes present in correlation matrix
        genes = gbmgenes.melt().dropna()
        genes = genes[genes.value.isin(corr_df1.index) & genes.value.isin(corr_df1.columns)]
        # select the columns and rows of df in gbmgenes
        corr_df = corr_df1.loc[genes.value,genes.value]
        # make a colour map for the genes
        lut = { col: f'tab:{clr}' for col, clr in zip(gbmgenes.columns,['red','blue'])}
        # map the colours to the genes
        row_colors = corr_df.index.map(genes.drop_duplicates(subset=['value']).set_index('value')['variable'])
        row_colors = row_colors.map(lut)
        row_colors = row_colors.fillna('white')
        # plot the heatmap of the correlation matrix
        g = sns.clustermap(data=corr_df, cmap='coolwarm',vmax=1, vmin=-1, row_colors=row_colors, col_colors=row_colors, row_cluster=False, col_cluster=False, dendrogram_ratio=0.11)
        # create a list of patches for the legend
        patches = [
            mpatches.Patch(color=color, label=column) for column, color in lut.items()
        ]
        # add legend for row_colours given by lut
        plt.legend(handles = patches , bbox_to_anchor=(1.7, 1), loc='upper left')
        g.fig.suptitle(suff, fontsize=30) 
        # save the figure
        plt.savefig('../figures/'+GSE+'/Correlation/'+GSM+'_Corrplot_'+suff+'.png',dpi=600)
        plt.close()
        plt.clf()

def wrapper(GSE,GSM):
    print(GSM)
    # Read the correlation matrix
    corr_df = pd.read_csv('../Output/'+GSE+'/Correlation/'+GSM+'_correlation.tsv',sep='\t',index_col=0)
    # Read the GBM genes
    gbmgenes = pd.read_csv("../Signatures/GBM_signatures.csv") 
    for suff in ['Nef','Ver']:
        gbmgenes1 = gbmgenes.filter(like=suff)
        corr_pair_plot(corr_df,gbmgenes1,GSE,GSM,suff)

# List of GSM IDs
files = glob.glob('../Output/'+GSE+'/Correlation/*_correlation.tsv')
GSMs = [os.path.basename(x).replace(f'_correlation.tsv','') for x in files]
# Make directory for figures
os.makedirs('../figures/'+GSE+'/Correlation/', exist_ok=True)

# Iterate over GSM samples
for GSM in GSMs:
    wrapper(GSE, GSM)