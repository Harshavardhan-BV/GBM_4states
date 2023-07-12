import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import itertools as it
plt.rcParams["svg.hashsalt"]=''

def corr_plot(corr_df, gbmgenes, GSE, GSM, suff):
    # Select only the columns that start with suff
    gbmgenes = gbmgenes.filter(like=suff)
    # Select only the genes present in correlation matrix
    genes = gbmgenes.melt().dropna()
    genes = genes[genes.value.isin(corr_df.index) & genes.value.isin(corr_df.columns)]
    # select the columns and rows of df in gbmgenes
    corr_df = corr_df.loc[genes.value,genes.value]
    # make a colour map for the genes
    lut = { col: f'tab:{clr}' for col, clr in zip(gbmgenes.columns,['red','orange','green','blue'])}
    # map the colours to the genes
    row_colors = corr_df.index.map(genes.drop_duplicates(subset=['value']).set_index('value')['variable'])
    row_colors = row_colors.map(lut)
    row_colors = row_colors.fillna('white')
    # plot the clustermap of the correlation matrix
    sns.clustermap(data=corr_df, cmap='coolwarm',vmax=1, vmin=-1, row_colors=row_colors, col_colors=row_colors)
    # add legend for row_colours given by lut
    # create a list of patches for the legend
    patches = [
        mpatches.Patch(color=color, label=column) for column, color in lut.items()
    ]
    plt.legend(handles = patches , bbox_to_anchor=(1.7, 1), loc='upper left')
    # save the figure
    plt.savefig('./figures/'+GSE+'/Correlation/'+GSM+'_Corrplot_'+suff+'.png',dpi=600)
    plt.close()
    plt.clf()

def corr_consistency(corr_df, gbmgenes, GSE, GSM, suff):
    # Select only the columns that start with suff
    gbmgenes = gbmgenes.filter(like=suff)
    combinats = list(it.combinations(gbmgenes.columns,2))
    cons = np.empty(len(combinats))
    for i in range(len(combinats)):
        combi = combinats[i]
        genes = np.empty(2,dtype='object')
        self_sum = np.empty(2)
        for j in range(2):
            gs = gbmgenes[combi[j]].dropna()
            # select only the genes that are in the correlation matrix
            gs = gs[gs.isin(corr_df.index) & gs.isin(corr_df.columns)]
            genes[j] = gs
            # subset the within group and across group correlations
            self = corr_df.loc[gs,gs]
            # calculate the mean within group correlation
            self_sum[j] = (np.tril(self).sum() - len(self)/2)/self.size
        cross = corr_df.loc[genes[0],genes[1]]
        # calculate the mean across group correlation
        cross_sum = cross.values.sum()/cross.size
        # Score the consistency
        cons[i] = (np.sum(self_sum) - cross_sum)/2
    df = pd.DataFrame(combinats,columns=['G1','G2'])
    df['Consistency'] = cons
    df.to_csv('./Output/'+GSE+'/Correlation/'+GSM+'_consistency_'+suff+'.tsv',sep='\t',index=False)
    plt.bar(range(len(combinats)),cons)
    plt.xticks(range(len(combinats)), combinats, rotation=45)
    plt.ylabel('Consistency')
    plt.tight_layout()
    plt.savefig('./figures/'+GSE+'/Correlation/'+GSM+'_Consistency_'+suff+'.svg')
    plt.clf()
    plt.close()
    
def wrapper(GSE,GSM):
    print(GSM)
    # Read the correlation matrix
    corr_df = pd.read_csv('./Output/'+GSE+'/Correlation/'+GSM+'_correlation.tsv',sep='\t',index_col=0)
    # Remove the NaNs
    corr_df.dropna(axis=0, inplace=True)
    # Read the GBM genes
    gbmgenes = pd.read_csv("./Signatures/GBM_signatures.csv") 
    for suff in ['Nef','Ver']:
        corr_plot(corr_df,gbmgenes,GSE,GSM,suff)
        corr_consistency(corr_df,gbmgenes,GSE,GSM,suff)

# GSE ID     
GSE = "GSE168004"
# GSE = "GSE131928"
# GSE = "GSE182109"
# GSE = "CCLE"

# List of GSM IDs
files = glob.glob('Output/'+GSE+'/Correlation/*_correlation.tsv')
GSMs = [os.path.basename(x).replace(f'_correlation.tsv','') for x in files]
# Make directory for figures
os.makedirs('figures/'+GSE+'/Correlation/', exist_ok=True)

# Iterate over GSM samples
for GSM in GSMs:
    wrapper(GSE, GSM)
