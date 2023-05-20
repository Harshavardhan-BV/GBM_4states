#%%
import numpy as np
import pandas as pd
#%%
df = pd.read_csv('Data/GSE182109_metadata.txt', sep = '\t')
#%%
df.loc[df.loc[:,"SubAssignment"]==np.nan, "SubAssignment"] = df.loc[df.loc[:,"SubAssignment"]==np.nan,"Assignment"]
## %%
# Get the frequency of each subassignment in each patient
nums = df.groupby('Patient')['Assignment'].value_counts()
# %%
nums.to_csv('./GSE182109_nums_Assignment.csv')
# %%
# Get the frequency of each subassignment in each patient
nums = df.groupby('Patient')[['Assignment','SubAssignment']].value_counts()
# %%
nums.to_csv('./GSE182109_nums_SubAssignment.csv')
# %%
