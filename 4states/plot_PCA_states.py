import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''

def loada_plot(GSE, GSM, suff):
    print(GSM,suff)
    load_df = pd.read_csv('./Output/'+GSE+'/PCA/'+GSM+'_loadings_'+suff+'.tsv',sep='\t',index_col=0)
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
    sns.clustermap(data=load_df, cmap='coolwarm',center=0, row_colors=row_colors, col_cluster=False)
    # create a list of patches for the legend
    patches = [
        mpatches.Patch(color=color, label=column) for column, color in lut.items()
    ]
    # add legend for row_colours given by lut
    plt.legend(handles = patches , bbox_to_anchor=(2, 1), loc='upper left')
    plt.savefig('./figures/'+GSE+'/PCA/'+GSM+'_loadings_plot_'+suff+'.svg')
    plt.close()
    plt.clf()

# GSE ID
GSE = "GSE168004"
# GSE = "GSE131928"
# GSE = "GSE182109"
# GSE = "CCLE"
# GSE = "TCGA"
# List files in Output
files = os.listdir('Output/'+GSE+'/PCA')
# Get the GSMs
GSMs = np.unique([x.rsplit('_',1)[0] for x in files])
os.makedirs('figures/'+GSE+'/PCA/', exist_ok=True)
# Iterative over the signatures
for suff in ['Nef','Ver']:
    # List files in Output
    files = glob.glob('Output/'+GSE+'/PCA/*_loadings_'+suff+'.tsv')
    # Get the GSMs
    GSMs = [os.path.basename(x).replace(f'_loadings_'+suff+'.tsv','') for x in files]
    # Iterate over GSM samples
    for GSM in GSMs:
        loada_plot(GSE,GSM, suff)
# %%
