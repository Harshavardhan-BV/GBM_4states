import os
import glob
import itertools as it
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''

def kde_cons(df, suff, maxnorm=False):
    if maxnorm:
        df['Consistency'] = df['Consistency']/df['Consistency'].max()
        suff+='_maxnorm'
    # Get combination from G1 and G2
    df['Combination'] = df['G1']+'-'+df['G2']
    # Plot the KDE
    sns.kdeplot(data=df, x='Consistency', hue='Combination', cut=0, fill=True, palette='tab10')
    plt.tight_layout()
    # Save the figure
    plt.savefig('./figures/Correlation/kdeplot_'+suff+'.svg')
    plt.clf()
    plt.close()

def box_cons(df, suff, maxnorm=False):
    if maxnorm:
        df['Consistency'] = df['Consistency']/df['Consistency'].max()
        suff+='_maxnorm'
    # Get combination from G1 and G2
    df['Combination'] = df['G1']+'-'+df['G2']
    # Plot the KDE
    sns.boxplot(data=df, x='Consistency', y='Combination')
    plt.tight_layout()
    # Save the figure
    plt.savefig('./figures/Correlation/boxplot_'+suff+'.svg')
    plt.clf()
    plt.close()

def corr_scraper(GSEs, suff):
    df_mega = pd.DataFrame()
    # Iterate through GSEs
    for GSE in GSEs:
        # List GSMs in each GSE
        files = glob.glob('./Output/'+GSE+'/Correlation/*_consistency_'+suff+'.tsv')
        GSMs = [os.path.basename(file).replace('_consistency_'+suff+'.tsv','') for file in files]
        # Iterate through GSMs
        for GSM in GSMs:
            # Read the correlation 
            df = pd.read_csv('./Output/'+GSE+'/Correlation/'+GSM+'_consistency_'+suff+'.tsv', sep='\t')
            # Add GSM and GSE columns
            df['GSM'] = GSM
            df['GSE'] = GSE
            # Append to mega dataframe
            df_mega = pd.concat([df_mega, df])
    return df_mega

os.makedirs('./figures/Correlation', exist_ok=True)

# Bulk datasets
datasets = pd.read_csv('./Datasets_Bulk.csv')
GSEs = datasets.GSE
# Neftel signatures
df = corr_scraper(GSEs, 'Nef')
box_cons(df, 'Bulk_Nef')
box_cons(df, 'Bulk_Nef', maxnorm=True)
# Verhaak signatures
df = corr_scraper(GSEs, 'Ver')
box_cons(df, 'Bulk_Ver')
box_cons(df, 'Bulk_Ver', maxnorm=True)
# Mix 
df = corr_scraper(GSEs, 'Mix')
box_cons(df, 'Bulk_Mix')
box_cons(df, 'Bulk_Mix', maxnorm=True)

# Single-cell datasets
datasets = pd.read_csv('./Datasets_SC.csv')
GSEs = datasets.GSE
# Neftel signatures
df = corr_scraper(GSEs, 'Nef')
box_cons(df, 'SC_Nef')
box_cons(df, 'SC_Nef', maxnorm=True)
# Verhaak signatures
box_cons(df, 'SC_Ver')
box_cons(df, 'SC_Ver', maxnorm=True)
# Mix 
box_cons(df, 'SC_Mix')
box_cons(df, 'SC_Mix', maxnorm=True)