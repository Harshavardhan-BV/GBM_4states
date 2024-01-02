#%%
import os
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''
sns.set_context('paper', font_scale=1.8)
#%%
# Read the GSE ID from command line
parser = argparse.ArgumentParser(
                    prog='hmap_Expr.py',
                    description='Plots heatmaps of the correlation between expression of genes in signatures')
parser.add_argument('GSE', type=str, nargs=1, help='GSE ID')
GSE = parser.parse_args().GSE[0]
#%%
def top_tfs(GSE, GSM, sig):
    # List the files
    # GSM3828672_NefMES-topTFs
    files = glob.glob('../Output/'+GSE+'/GRN/'+GSM+'_'+sig+'*-topTFs.tsv')
    top = []
    # Iterate through the files
    for file in files:
        # Get the signature name
        score = os.path.basename(file).replace(GSM+'_','').replace('-topTFs.tsv','')
        # Read the top TFs and add to top
        top += pd.read_csv(file, sep='\t', header=None).values.flatten().tolist()
    return top
#%%
def corr_plot(GSE, GSM, suff):
    corr_df = pd.read_csv('../Output/'+GSE+'/GRN/'+GSM+'_TF-corr.tsv', sep='\t', index_col=0)
    # Read the top TFs from the files
    files = glob.glob('../Output/'+GSE+'/GRN/'+GSM+'_'+suff+'*-topTFs.tsv')
    top = []
    row_colors = []
    palette = {'NefNPC':'tab:blue', 'NefOPC':'tab:green', 'NefAC':'tab:orange', 'NefMES':'tab:red', 'VerPN':'tab:blue', 'VerNL':'tab:green', 'VerCL':'tab:orange', 'VerMES':'tab:red'}
    # Iterate through the files
    for file in files:
        # Get the signature name
        score = os.path.basename(file).replace(GSM+'_','').replace('-topTFs.tsv','')
        # Read the top TFs and add to top
        top_df = pd.read_csv(file, sep='\t', header=None).values.flatten().tolist()
        row_colors += len(top_df)*[palette[score]]
        top += top_df
    # Select the top TFs from the correlation dataframe
    corr_df = corr_df.loc[top,top]
    # plot the clustermap of the correlation matrix
    sns.clustermap(data=corr_df, cmap='coolwarm',vmax=1, vmin=-1,row_colors=row_colors,col_colors=row_colors)
    # save the figure
    plt.savefig('../figures/'+GSE+'/GRN/'+GSM+suff+'_TF-corrplot.svg')
    plt.close()
    plt.clf()
#%%
# List of GSM IDs
files = glob.glob('../Output/'+GSE+'/GRN/*_TF-corr.tsv')
GSMs = [os.path.basename(x).replace(f'_TF-corr.tsv','') for x in files]
# Make directory for figures
os.makedirs('../figures/'+GSE+'/GRN/', exist_ok=True)
#%%
# Iterate over GSM samples
for GSM in GSMs:
    corr_plot(GSE, GSM, suff='Nef')
    corr_plot(GSE, GSM, suff='Ver')
# %%
