import os
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

# Read the GSE ID and sc from command line
parser = argparse.ArgumentParser(
                    prog='Corr_split_GSEA.py',
                    description='Calculates correlation and pvalues between the AUCell/ssGSEA scores for a given GSE')
parser.add_argument('GSE', type=str, nargs=1, help='GSE ID')
parser.add_argument('-sc', action=argparse.BooleanOptionalAction, help='Add this flag if single cell data')
GSE = parser.parse_args().GSE[0]
sc = parser.parse_args().sc

def spearman_pval(x,y):
    return spearmanr(x,y)[1]

def hmap_corr(df,ax, title):
    corr_df = df.corr(method='spearman')
    pval_df = df.corr(method=spearman_pval)
    np.fill_diagonal(pval_df.values, 0)
    pval_df = pval_df.applymap(lambda x: '*' if x < 0.05 else '')
    sns.heatmap(corr_df, cmap='coolwarm', vmax=1, vmin=-1, annot=pval_df, fmt='', annot_kws={"size": 20}, ax=ax, cbar=False)
    ax.set_title(title+' n='+str(len(df)))

def aucorr(GSE, GSM, suff):
    print(GSM)
    df = pd.read_csv('../Output/'+GSE+'/'+score+'/'+GSM+'-'+score+'.csv', index_col=0)
    if suff=='Nef':
        extm = ['NefMES','NefNPC']
        intr = ['NefAC','NefOPC']
    elif suff=='Ver':
        extm = ['VerMES','VerPN']
        intr = ['VerCL','VerNL']
    # Select only the Ver or Nef signatures
    df = df.loc[:,extm+intr]
    # Plot the heatmap
    fig , ax = plt.subplots(1,3, figsize=(30,10))
    # Plot all the correlations
    hmap_corr(df,ax[0], 'All')
    # Split the dataframe into extreme and intermediate cases
    df_extm = df[df.idxmax(axis=1).isin(extm)]
    df_intr = df[df.idxmax(axis=1).isin(intr)]
    # Plot the heatmap for extreme and intermediate cases
    hmap_corr(df_extm,ax[1], 'Extreme')
    hmap_corr(df_intr,ax[2], 'Intermediate')
    # Save the figure
    plt.tight_layout()
    plt.savefig('../figures/'+GSE+'/GSEA/GSEA_'+GSM+'_splitcorr_'+suff+'.svg')
    plt.clf()
    plt.close()

# List of GSM IDs
score = 'AUCell' if sc else 'ssGSEA'
files = glob.glob('../Output/'+GSE+'/'+score+'/*-'+score+'.csv')
GSMs = [os.path.basename(x).replace(f'-'+score+'.csv','') for x in files]
# Make directory for output
os.makedirs('../figures/'+GSE+'/GSEA/', exist_ok=True)

# Iterate over GSM samples
for GSM in GSMs:
    aucorr(GSE,GSM,'Nef')
    aucorr(GSE,GSM,'Ver')