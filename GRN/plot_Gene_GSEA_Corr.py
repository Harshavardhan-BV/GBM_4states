import os
import glob
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams["svg.hashsalt"]=''
sns.set_context('paper', font_scale=1.8)

# Read the GSE ID from command line
parser = argparse.ArgumentParser(
                    prog='bleh.py',
                    description='blehhh')
parser.add_argument('GSE', type=str, nargs=1, help='GSE ID')
parser.add_argument('-R', '--R_cutoff', type=float, nargs=1, help='Cutoff of correlation coefficient', default=[0.3])
GSE = parser.parse_args().GSE[0]
R_cut = parser.parse_args().R_cutoff[0]

def top_TF(GSE, GSM, x, top=10, save=True):
    # Read the correlation dataframe
    corr_df = pd.read_csv('../Output/'+GSE+'/GRN/'+GSM+'_GG-corr.tsv', sep='\t', index_col=0)
    # Read the list of transcription factors
    tfs = pd.read_csv('../Signatures/Human_transcription_factors.txt', sep='\t')
    # Get the intersection of the two dataframes
    corr_df = corr_df[corr_df.index.isin(tfs['tfs'])]
    # Get the top genes
    top_df = corr_df.sort_values(by=x, ascending=False).iloc[:top,:].index.to_frame()
    # Save the top genes to a file
    top_df.to_csv('../Output/'+GSE+'/GRN/'+GSM+'_'+x+'-topTFs.tsv', sep='\t', index=False, header=False)
    return top_df.values.flatten().tolist()
    
def corr_vs_plot(GSE, GSM, x, y, tf=False, top=None):
    # Read the correlation dataframe
    corr_df = pd.read_csv('../Output/'+GSE+'/GRN/'+GSM+'_GG-corr.tsv', sep='\t', index_col=0)
    if tf:
        # Read the list of transcription factors
        tfs = pd.read_csv('../Signatures/Human_transcription_factors.txt', sep='\t')
        # Get the intersection of the two dataframes
        corr_df = corr_df[corr_df.index.isin(tfs['tfs'])]
        fname = '../figures/'+GSE+'/GRN/'+GSM+'_'+x+'-'+y+'_GTF-corr.svg'
    else:
        fname = '../figures/'+GSE+'/GRN/'+GSM+'_'+x+'-'+y+'_GG-corr.svg'
    # plot the scatter of x vs y
    plt.figure(figsize=(6,6))
    sns.scatterplot(data=corr_df, x=x, y=y, c='k', s=10, edgecolor=None)
    # Add a horizontal line at y=0
    plt.axhline(y=0, color='r', linestyle='--')
    # Add a vertical line at x=0
    plt.axvline(x=0, color='r', linestyle='--')
    # Add a grey box around 0.2 and -0.2
    plt.axhspan(-0.2, 0.2, color='grey', alpha=0.2)
    plt.axvspan(-0.2, 0.2, color='grey', alpha=0.2)
    # Annotate the genes given in the genelist if top is not None
    if top is not None:
        from adjustText import adjust_text
        # Identify the index of top genes
        top_df = corr_df.loc[top,:]
        # Annotate the top genes
        texts = [plt.text(top_df.loc[i,x], top_df.loc[i,y], i, c='b', size='x-small') for i in top]
        adjust_text(texts, force_text=(0.4, 0.7), arrowprops=dict(arrowstyle="-", color='b', lw=0.5) )
    # save the figure
    plt.tight_layout()
    plt.savefig(fname)
    plt.close()
    plt.clf()

def TF_cluster(GSE, GSM, R_cut=0.3):
    # Read the correlation dataframe
    corr_df = pd.read_csv('../Output/'+GSE+'/GRN/'+GSM+'_GG-corr.tsv', sep='\t', index_col=0)
    # Read the list of transcription factors
    tfs = pd.read_csv('../Signatures/Human_transcription_factors.txt', sep='\t')
    # Get the intersection of the two dataframes
    corr_df = corr_df[corr_df.index.isin(tfs['tfs'])]
    # Select values where any of the correlation coefficient is greater than R_cut
    corr_df = corr_df[(corr_df.abs() > R_cut).any(axis=1)]
    # plot the clustermap of the correlation matrix
    sns.clustermap(data=corr_df, cmap='coolwarm',center=0, figsize=(20,20))
    # save the figure
    plt.savefig('../figures/'+GSE+'/GRN/'+GSM+'_GTF-cluster.svg')
    plt.close()
    plt.clf()

# List of GSM IDs
files = glob.glob('../Output/'+GSE+'/GRN/*_GG-corr.tsv')
GSMs = [os.path.basename(x).replace(f'_GG-corr.tsv','') for x in files]
# Make directory for figures
os.makedirs('../figures/'+GSE+'/GRN/', exist_ok=True)

# Iterate over GSM samples
for GSM in GSMs:
    print(GSM)
    # TF_cluster(GSE, GSM, R_cut)
    x = [['NefNPC', 'NefMES'],['VerPN', 'VerMES']]
    for xi in x:
        corr_vs_plot(GSE, GSM, xi[0], xi[1])
        topgenes = []
        for xj in xi:
            topgenes+=top_TF(GSE, GSM, xj, top=10)
        corr_vs_plot(GSE, GSM,  xi[0], xi[1], tf=True, top=topgenes)