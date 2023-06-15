#%%
import numpy as np
import pandas as pd
#%%
df = pd.read_csv('Data/GSE182109_metadata.txt', sep = '\t')
# %%
mapper ={   'h-microglia' : 'Mic', 
            'i-microglia':'Mic',
            'a-microglia':'Mic',
            's-mac 1':'Mac',
            's-mac 2':'Mac',
            'CD8 TCells' : 'CD8T'
        }
#%%
# Replace the values in the SubAssignment column with the index in the mapper
df['NewAssignment'] = df['SubAssignment'].map(mapper)
# %%
df.loc[df['Assignment']=='Glioma','NewAssignment'] = 'Glioma'
# %%
df.to_csv('Data/GSE182109_metadata.tsv', sep = '\t', )
# %%
