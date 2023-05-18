#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#%%
df = pd.read_csv('./Output/CCLE/ssgsea_GBM/gseapy.gene_set.ssgsea.report.csv')
# %%
df = df.pivot(index='Name',columns='Term',values='NES')
#%%
sns.pairplot(df)
plt.savefig('./figures/CCLE/ssgsea_pairplot_all.svg')
# %%
mesdf = df[['MES','MES1','MES2']]
# %%
npcdf = df[['NPC','NPC1','NPC2']]
# %%
fourdf = df[['NPC','AC','OPC','MES']]
# %%
sns.pairplot(mesdf)
plt.savefig('./figures/CCLE/ssgsea_pairplot_mes.svg')
#%%
sns.pairplot(npcdf)
plt.savefig('./figures/CCLE/ssgsea_pairplot_npc.svg')
# %%
sns.pairplot(fourdf)
plt.savefig('./figures/CCLE/ssgsea_pairplot_4D.svg')
# %%
sns.clustermap(fourdf.corr(), cmap='coolwarm',vmax=1, vmin=-1)
plt.savefig('./figures/CCLE/ssgsea_corr_4D.svg')
# %%
