import os
import glob
import itertools as it
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''

def kde_cons(df, suff):
    # Get combination from G1 and G2
    df['Combination'] = df['G1']+'-'+df['G2']
    # Plot the KDE
    sns.kdeplot(data=df, x='Consistency', hue='Combination', cut=0, fill=True, palette='tab10')
    plt.tight_layout()
    # Save the figure
    plt.savefig('./figures/Correlation/kdeplot_'+suff+'.svg')
    plt.clf()
    plt.close()

def box_cons(df, suff):
    # Get combination from G1 and G2
    df['Combination'] = df['G1']+'-'+df['G2']
    # Plot the KDE
    sns.boxplot(data=df, x='Consistency', y='Combination')
    plt.tight_layout()
    # Save the figure
    plt.savefig('./figures/Correlation/boxplot_'+suff+'.svg')
    plt.clf()
    plt.close()

def corr_scraper(GSEs, suff, maxnorm=False):
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
            if maxnorm:
                df['Consistency'] = df['Consistency'] / df['Consistency'].max()
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
df = corr_scraper(GSEs, 'Nef', maxnorm=True)
box_cons(df, 'Bulk_Nef_maxnorm')
# Verhaak signatures
df = corr_scraper(GSEs, 'Ver')
box_cons(df, 'Bulk_Ver')
df = corr_scraper(GSEs, 'Nef', maxnorm=True)
box_cons(df, 'Bulk_Ver_maxnorm')
# Mix 
df = corr_scraper(GSEs, 'Mix')
box_cons(df, 'Bulk_Mix')
df = corr_scraper(GSEs, 'Mix', maxnorm=True)
box_cons(df, 'Bulk_Mix_maxnorm')


# Single-cell datasets
datasets = pd.read_csv('./Datasets_SC.csv')
GSEs = datasets.GSE
# Neftel signatures
df = corr_scraper(GSEs, 'Nef')
box_cons(df, 'SC_Nef')
df = corr_scraper(GSEs, 'Nef', maxnorm=True)
box_cons(df, 'SC_Nef_maxnorm')
# Verhaak signatures
df = corr_scraper(GSEs, 'Ver')
box_cons(df, 'SC_Ver')
df = corr_scraper(GSEs, 'Nef', maxnorm=True)
box_cons(df, 'SC_Ver_maxnorm')
# Mix 
df = corr_scraper(GSEs, 'Mix')
box_cons(df, 'SC_Mix')
df = corr_scraper(GSEs, 'Mix', maxnorm=True)
box_cons(df, 'SC_Mix_maxnorm')