#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#%%

def corr_plot(GSE, GSM):
    corr_df = pd.read_csv('./Output/'+GSE+'/'+GSM+'_correlation.tsv',sep='\t',index_col=0)
    corr_df.isna()
    corr_df.dropna(axis=0, inplace=True)
    # %%
    sns.clustermap(corr_df, cmap='coolwarm',vmax=1, vmin=-1)
    plt.savefig('./figures/'+GSE+'/'+GSM+'_Corrplot.svg')
    #%%
    gbmgenes = pd.read_excel("./Signatures/1-s2.0-S0092867419306877-mmc2.xlsx",skiprows=4)
    gbmgenes = gbmgenes.iloc[:,0:6]
    # %%
    # make a colour map for the genes
    lut = {
        'MES1': 'tab:red',
        'MES2': 'tab:red',
        'NPC1': 'tab:blue',
        'NPC2': 'tab:blue',
        'OPC': 'tab:green',
        'AC': 'tab:orange',
    }
    row_colors = corr_df.index.map(gbmgenes.melt().drop_duplicates(subset=['value']).set_index('value')['variable'])
    #%%
    row_colors = row_colors.map(lut)
    row_colors = row_colors.fillna('white')
    #%%
    sns.clustermap(data=corr_df, cmap='coolwarm',vmax=1, vmin=-1, row_colors=row_colors)
    plt.savefig('./figures/'+GSE+'/'+GSM+'_Corrplot_coloured.svg')
    plt.close()
    plt.clf()
# %%
GSE = "GSE131928"
GSMs = ["GSM3828672","GSM3828673"]
#%%
for GSM in GSMs:
    corr_plot(GSE, GSM)
# %%
