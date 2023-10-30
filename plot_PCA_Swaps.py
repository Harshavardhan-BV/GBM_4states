#%%
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''
#%%

def pc1_swapbox(GSE, GSM, suff):
    df = pd.read_csv('Output/'+GSE+'/PCA_Swaps/'+GSM+'_'+suff+'_swaps.csv')
    df['PC1_Var'] = df['PC1_Var']*100
    sns.boxplot(x='Swap_pct', y='PC1_Var',hue='Signature', data=df)
    plt.xlabel('Percentage of swapped genes')
    plt.ylabel('Variance explained by PC1')
    plt.savefig('figures/'+GSE+'/PCA_Swaps/'+GSM+'_'+suff+'_PC1_Swaps.svg')
    plt.close()
    plt.clf()

GSE = "GSE131928"
# GSE = "TCGA"

# Make the output directory
os.makedirs('figures/'+GSE+'/PCA_Swaps/', exist_ok=True)
# Iterative over the signatures
for suff in ['Nef','Ver']:
    # List files in Output
    files = glob.glob('Output/'+GSE+'/PCA_Swaps/*_'+suff+'_swaps.csv')
    # Get the GSMs
    GSMs = [os.path.basename(x).replace(f'_'+suff+'_swaps.csv','') for x in files]
    # Iterate over GSM samples
    for GSM in GSMs:
        pc1_swapbox(GSE,GSM, suff)
