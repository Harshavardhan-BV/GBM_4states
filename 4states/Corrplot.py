import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''

def corr_plot(GSE, GSM):
    # Read the correlation matrix
    corr_df = pd.read_csv('./Output/'+GSE+'/Correlation/'+GSM+'_correlation.tsv',sep='\t',index_col=0)
    # Remove the NaNs
    corr_df.dropna(axis=0, inplace=True)
    # Read the GBM genes
    gbmgenes = pd.read_excel("./Signatures/1-s2.0-S0092867419306877-mmc2.xlsx",skiprows=4)
    gbmgenes = gbmgenes.iloc[:,0:6]
    # make a colour map for the genes
    lut = {
        'MES1': 'tab:red',
        'MES2': 'tab:red',
        'NPC1': 'tab:blue',
        'NPC2': 'tab:blue',
        'OPC': 'tab:green',
        'AC': 'tab:orange',
    }
    # map the colours to the genes
    row_colors = corr_df.index.map(gbmgenes.melt().drop_duplicates(subset=['value']).set_index('value')['variable'])
    row_colors = row_colors.map(lut)
    row_colors = row_colors.fillna('white')
    # plot the clustermap of the correlation matrix
    sns.clustermap(data=corr_df, cmap='coolwarm',vmax=1, vmin=-1, row_colors=row_colors)
    # save the figure
    plt.savefig('./figures/'+GSE+'/'+GSM+'_Corrplot_coloured.svg')
    plt.close()
    plt.clf()

# GSE ID     
GSE = "GSE168004"
# List of GSM IDs
GSMs = [x.replace("_correlation.tsv",'') for x in os.listdir('Output/'+GSE+'/Correlation/')]
# Make directory for figures
os.makedirs('figures/'+GSE+'/Correlation/', exist_ok=True)

# Iterate over GSM samples
for GSM in GSMs:
    corr_plot(GSE, GSM)
