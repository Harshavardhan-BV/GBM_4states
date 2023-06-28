#%%
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#%%

def loada_plot(GSE, GSM):
    print(GSM)
    load_df = pd.read_csv('./Output/'+GSE+'/PCA/'+GSM+'_loadings.tsv',sep='\t',index_col=0)
    gbmgenes = pd.read_excel("./Signatures/1-s2.0-S0092867419306877-mmc2.xlsx",skiprows=4)
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
    plt.savefig('./figures/'+GSE+'/PCA/'+GSM+'_loadings_plot.svg')
    plt.close()
    plt.clf()
# %%
# GSE ID
GSE = 'GSE168004'
# List files in Output
files = os.listdir('Output/'+GSE+'/PCA')
# Get the GSMs
GSMs = np.unique([x.rsplit('_',1)[0] for x in files])
os.makedirs('figures/'+GSE+'/PCA/', exist_ok=True)
#%%
for GSM in GSMs:
    loada_plot(GSE, GSM)
# %%
