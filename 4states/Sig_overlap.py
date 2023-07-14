#%%
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''
# Read the signatures
df = pd.read_csv('./Signatures/GBM_signatures.csv')
# Make directory for figures and output
os.makedirs('./figures/Sig_overlap/', exist_ok=True)
os.makedirs('./Output/Sig_overlap/', exist_ok=True)
# Create a nxn matrix for overlap
n = len(df.columns)
overlap = np.zeros((n,n))
# Iterate over the signatures
for i in range(n):
    # Select the signature
    i_col = df.iloc[:,i].dropna()
    # Overlap with itself to be number of genes in the signature
    overlap[i,i] = len(i_col)
    for j in range(i+1, n):
        # Select the signature
        j_col = df.iloc[:,j].dropna()
        # Intersect the two signatures
        inter = np.intersect1d(i_col, j_col)
        # Save the overlapping genes if any
        if len(inter)>0:
            np.savetxt('./Output/Sig_overlap/'+df.columns[i]+'_'+df.columns[j]+'.txt', inter, fmt='%s')
        # Store the number of genes in overlap matrix
        overlap[i,j] = len(inter)
# Making it symmetric
overlap = overlap+overlap.T-np.diag(np.diag(overlap))
# Convert to dataframe
overlap = pd.DataFrame(overlap, index=df.columns, columns=df.columns)

def overlap_plt(overlap, suff):
    # Filter the overlap for signature set
    if suff!='all':
        overlap = overlap.filter(like=suff, axis=0).filter(like=suff, axis=1)
    # Plot a heatmap
    sns.heatmap(overlap, annot=True, cmap='viridis',fmt=".0f")
    plt.tight_layout()
    # Save the figure
    plt.savefig('./figures/Sig_overlap/signature_overlap_'+suff+'.svg')
    plt.clf()
    plt.close()

overlap_plt(overlap, 'all')
overlap_plt(overlap, 'Nef')
overlap_plt(overlap, 'Ver')