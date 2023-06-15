#%%
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams["svg.hashsalt"]=''
#%%
def corr_plot_AUCell(GSE,GSM):
    print(GSM)
    df = pd.read_csv('Output/'+GSE+'/'+GSM+'-AUCell.csv', index_col=0)
    df = df.loc[:, (df != 0).any(axis=0)]
    mesdf = df[['MES','MES1','MES2']]
    npcdf = df[['NPC','NPC1','NPC2']]
    fourdf = df[['NPC','AC','OPC','MES']]
    ns.clustermap(df.corr(), cmap='coolwarm',vmax=1, vmin=-1)
    plt.savefig('figures/'+GSE+'/AUCell_'+GSM+'_corrplot_all.svg')
    plt.clf()
    plt.close()
    ns.clustermap(mesdf.corr(), cmap='coolwarm',vmax=1, vmin=-1)
    plt.savefig('figures/'+GSE+'/AUCell_'+GSM+'_corrplot_mes.svg')
    plt.clf()
    plt.close()
    ns.clustermap(npcdf.corr(), cmap='coolwarm',vmax=1, vmin=-1)
    plt.savefig('figures/'+GSE+'/AUCell_'+GSM+'_corrplot_npc.svg')
    plt.clf()
    plt.close()
    ns.clustermap(fourdf.corr(), cmap='coolwarm',vmax=1, vmin=-1)
    plt.savefig('figures/'+GSE+'/AUCell_'+GSM+'_corrplot_4D.svg')
    plt.clf()
    plt.close()

#%%        
# GSE ID
GSE = "GSE168004"
GSMs = ["OSM_celllines", "mgg23", "mgg75"]
#%%
os.makedirs('figures/'+GSE+'/AUCell/', exist_ok=True)
for GSM in GSMs:
    corr_plot_AUCell(GSE,GSM)
