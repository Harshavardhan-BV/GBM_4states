import os
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''

# Read the GSE ID from command line
parser = argparse.ArgumentParser(
                    prog='varplot_PCA_Swaps.py',
                    description='Plots the variance explained of PC1 for different percentages of swapped genes')
parser.add_argument('GSE', type=str, nargs=1, help='GSE ID')
GSE = parser.parse_args().GSE[0]

def pc1_swapbox(GSE, GSM, suff):
    df = pd.read_csv('../Output/'+GSE+'/PCA_Swaps/'+GSM+'_'+suff+'_swaps.csv')
    df['PC1_Var'] = df['PC1_Var']*100
    sns.boxplot(x='Swap_pct', y='PC1_Var',hue='Signature', data=df)
    plt.xlabel('Percentage of swapped genes')
    plt.ylabel('Variance explained by PC1')
    plt.savefig('../figures/'+GSE+'/PCA_Swaps/'+GSM+'_'+suff+'_PC1_Swaps.svg')
    plt.close()
    plt.clf()

# Make the output directory
os.makedirs('../figures/'+GSE+'/PCA_Swaps/', exist_ok=True)
# Iterative over the signatures
for suff in ['Nef','Ver']:
    # List files in Output
    files = glob.glob('../Output/'+GSE+'/PCA_Swaps/*_'+suff+'_swaps.csv')
    # Get the GSMs
    GSMs = [os.path.basename(x).replace(f'_'+suff+'_swaps.csv','') for x in files]
    # Iterate over GSM samples
    for GSM in GSMs:
        pc1_swapbox(GSE,GSM, suff)
