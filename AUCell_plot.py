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
    sns.pairplot(df, plot_kws = dict(edgecolor="none"))
    plt.savefig('figures/'+GSE+'/AUCell_'+GSM+'_pairplot_all.svg')
    plt.clf()
    plt.close()
    sns.clustermap(df.corr(), cmap='coolwarm',vmax=1, vmin=-1)
    plt.savefig('figures/'+GSE+'/AUCell_'+GSM+'_corrplot_all.svg')
    plt.clf()
    plt.close()
    sns.pairplot(mesdf, plot_kws = dict(edgecolor="none"))
    plt.savefig('figures/'+GSE+'/AUCell_'+GSM+'_pairplot_mes.svg')
    plt.clf()
    plt.close()
    sns.clustermap(mesdf.corr(), cmap='coolwarm',vmax=1, vmin=-1)
    plt.savefig('figures/'+GSE+'/AUCell_'+GSM+'_corrplot_mes.svg')
    plt.clf()
    plt.close()
    sns.pairplot(npcdf, plot_kws = dict(edgecolor="none"))
    plt.savefig('figures/'+GSE+'/AUCell_'+GSM+'_pairplot_npc.svg')
    plt.clf()
    plt.close()
    sns.clustermap(npcdf.corr(), cmap='coolwarm',vmax=1, vmin=-1)
    plt.savefig('figures/'+GSE+'/AUCell_'+GSM+'_corrplot_npc.svg')
    plt.clf()
    plt.close()
    sns.pairplot(fourdf, plot_kws = dict(edgecolor="none"))
    plt.savefig('figures/'+GSE+'/AUCell_'+GSM+'_pairplot_4D.svg')
    plt.clf()
    plt.close()
    sns.clustermap(fourdf.corr(), cmap='coolwarm',vmax=1, vmin=-1)
    plt.savefig('figures/'+GSE+'/AUCell_'+GSM+'_corrplot_4D.svg')
    plt.clf()
    plt.close()

#%%        
# GSE ID
GSE = "GSE131928"
GSMs = ["GSM3828672","GSM3828673"]
#%%
os.makedirs('figures/'+GSE+'/AUCell/', exist_ok=True)
for GSM in GSMs:
    corr_plot_AUCell(GSE,GSM)