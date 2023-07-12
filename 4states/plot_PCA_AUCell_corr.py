import os
import glob
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams["svg.hashsalt"]=''

def pca_corr_plot(GSE,GSM, suff):
    print(GSM,suff)
    try:
        # Read the correlation files
        df_top = pd.read_csv('Output/'+GSE+'/PCA-full/'+GSM+'_top50-corr_'+suff+'.tsv', sep='\t', index_col=0)
        df_last = pd.read_csv('Output/'+GSE+'/PCA-full/'+GSM+'_last50-corr_'+suff+'.tsv', sep='\t', index_col=0)
    except FileNotFoundError:
        # If no correlation file, skip
        print('No correlation file for '+GSM+'found')
        return 
    # Make a heatmap for top and last seperately
    fig, ax = plt.subplots(1,2, figsize=(10,10))
    sns.heatmap(df_top, cmap='coolwarm',vmax=1,vmin=-1, ax = ax[0])
    sns.heatmap(df_last, cmap='coolwarm',vmax=1,vmin=-1, ax = ax[1])
    plt.tight_layout()
    # Save the figure
    plt.savefig('figures/'+GSE+'/PCA-full/'+GSM+'_load-corr_'+suff+'.svg')
    plt.clf()
    plt.close()
      
# GSE ID
GSE = "GSE168004"
# GSE = "GSE131928"
# GSE = "GSE182109"
# GSE = "CCLE"

# Make directory for figures
os.makedirs('figures/'+GSE+'/PCA-full/', exist_ok=True)
# Iterative over the signatures
for suff in ['Nef','Ver']:
    # List files in Output
    files = glob.glob('Output/'+GSE+'/PCA-full/*_top50-corr_'+suff+'.tsv')
    # Get the GSMs
    GSMs = [os.path.basename(x).replace(f'_top50-corr_'+suff+'.tsv','') for x in files]
    # Iterate over GSM samples
    for GSM in GSMs:
        pca_corr_plot(GSE,GSM, suff)
