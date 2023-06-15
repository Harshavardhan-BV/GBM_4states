#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#%%
df = pd.read_excel('Signatures/1-s2.0-S0092867419306877-mmc2.xlsx', skiprows=4)
# %%
df = df.iloc[:,0:6]
n = len(df.columns)
overlap = np.zeros((n,n))
# %%
for i in range(n):
    i_col = df.iloc[:,i].dropna()
    overlap[i,i] = len(i_col)
    for j in range(i+1, n):
        j_col = df.iloc[:,j].dropna()
        inter = np.intersect1d(i_col, j_col)
        overlap[i,j] = len(inter)
        
overlap = overlap+overlap.T-np.diag(np.diag(overlap))
# %%
overlap = pd.DataFrame(overlap, index=df.columns, columns=df.columns)
# %%
sns.heatmap(overlap, annot=True, cmap='viridis')
plt.savefig('./figures/signature_overlap.svg')
# %%
