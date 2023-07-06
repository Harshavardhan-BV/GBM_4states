import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams["svg.hashsalt"]=''

def corr_plot_AUCell(GSE,GSM):
    print(GSM)
    df = pd.read_csv('Output/'+GSE+'/AUCell/'+GSM+'-AUCell.csv', index_col=0)
    df = df.loc[:, (df != 0).any(axis=0)]
    fourdf = df[['NPC','AC','OPC','MES']]
    sns.clustermap(df.corr(method='spearman'), cmap='coolwarm',vmax=1, vmin=-1)
    plt.savefig('figures/'+GSE+'/AUCell/AUCell_'+GSM+'_corrplot_all.svg')
    plt.clf()
    plt.close()
    sns.clustermap(fourdf.corr(method='spearman'), cmap='coolwarm',vmax=1, vmin=-1)
    plt.savefig('figures/'+GSE+'/AUCell/AUCell_'+GSM+'_corrplot_4D.svg')
    plt.clf()
    plt.close()
     
# GSE ID
GSE = "GSE168004"
# List of GSM IDs
GSMs = [x.replace("-AUCell.csv",'') for x in os.listdir('Output/'+GSE+'/AUCell/')]
# Make directory for figures
os.makedirs('figures/'+GSE+'/AUCell/', exist_ok=True)

# Iterate over GSM samples
for GSM in GSMs:
    corr_plot_AUCell(GSE,GSM)
