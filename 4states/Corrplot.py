#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#%%
counts = pd.read_csv('./Data/CCLE/OmicsExpressionProteinCodingGenesTPMLogp1.csv',index_col=0)
#%%
celllines = pd.read_csv('./Data/CCLE/cell lines in Glioblastoma.csv')
celllines = celllines['Depmap Id'].tolist()
#%%
gbmgenes = pd.read_excel("./Signatures/1-s2.0-S0092867419306877-mmc2.xlsx",skiprows=4)
gbmgenes = gbmgenes.iloc[:,0:6]
# %%
# Get the values of counts for index in celllines
counts = counts[counts.index.isin(celllines)]
# %%
# removes the part within () for columns
counts.columns = counts.columns.str.replace(r" \(.*\)","",regex=True)
#%% Export the counts for the GBM celllines only
counts.T.to_csv('./Data_generated/CCLE/GBMCounts.tsv',sep='\t')
# %%
genes = gbmgenes.melt().dropna().value.tolist()
# %%
# Get the common elements between genes and counts.columns
genes = counts.columns.intersection(genes)
# %%
counts = counts.loc[:,genes]
# %%
corr_df = counts.corr()
# %%
sns.clustermap(corr_df, cmap='coolwarm',vmax=1, vmin=-1)
plt.savefig('./figures/CCLE/Corrplot.svg')
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

row_colors = counts.columns.map(gbmgenes.melt().drop_duplicates(subset=['value']).set_index('value')['variable'])
#%%
row_colors = row_colors.map(lut)
#%%
sns.clustermap(data=corr_df, cmap='coolwarm',vmax=1, vmin=-1, row_colors=row_colors)
plt.savefig('./figures/CCLE/Corrplot_coloured.svg')

# %%
