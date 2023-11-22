import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read the GSE ID from command line
parser = argparse.ArgumentParser(
                    prog='loadingsplot_PCA_Corr.py',
                    description='Plots a barplot of the loadings of the PCA on correlation matrix')
parser.add_argument('GSE', type=str, nargs=1, help='GSE ID')
GSE = parser.parse_args().GSE[0]

def loada_plot(GSE, GSM):
    print(GSM)
    load_df = pd.read_csv('../Output/'+GSE+'/Corr-PCA/'+GSM+'_loadings.tsv',sep='\t',index_col=0)
    gbmgenes = pd.read_excel("../Signatures/1-s2.0-S0092867419306877-mmc2.xlsx",skiprows=4)
    gbmgenes = gbmgenes.iloc[:,0:6]
    # select only the genes in the signature
    keys = gbmgenes.melt()['value'].unique()
    keys = keys[np.isin(keys,load_df.index)]
    load_df = load_df.loc[keys,:]
    # make a colour map for the genes
    lut = {
        'MES1': 'tab:red',
        'MES2': 'tab:red',
        'NPC1': 'tab:blue',
        'NPC2': 'tab:blue',
        'OPC': 'tab:green',
        'AC': 'tab:orange',
    }
    row_colors = load_df.index.map(gbmgenes.melt().drop_duplicates(subset=['value']).set_index('value')['variable'])
    row_colors = row_colors.map(lut)
    row_colors = row_colors.fillna('white')
    # Plot the loadings matrix
    sns.clustermap(data=load_df, cmap='coolwarm',center=0, row_colors=row_colors, col_cluster=False)
    plt.savefig('../figures/'+GSE+'/Corr-PCA/'+GSM+'_loadings_plot.svg')
    plt.close()
    plt.clf()

# List files in Output
files = os.listdir('../Output/'+GSE+'/Corr-PCA')
# Get the GSMs
GSMs = np.unique([x.rsplit('_',1)[0] for x in files])
os.makedirs('../figures/'+GSE+'/Corr-PCA/', exist_ok=True)

for GSM in GSMs:
    loada_plot(GSE, GSM)

