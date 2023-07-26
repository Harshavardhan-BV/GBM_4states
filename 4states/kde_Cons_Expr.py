import os
import glob
import itertools as it
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''

GSEs = ['CCLE','TCGA','GSE131928', 'GSE168004','GSE182109']

def kde_cons(GSEs, suff):
    # Read the consistency across all GSEs
    df = corr_scraper(GSEs, suff)
    # Get combination from G1 and G2
    df['Combination'] = df['G1']+'-'+df['G2']
    # Plot the KDE
    sns.kdeplot(data=df, x='Consistency', hue='Combination', cut=0, fill=True, palette='tab10')
    plt.tight_layout()
    # Save the figure
    plt.savefig('./figures/Correlation/kdeplot_'+suff+'.svg')
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
# Neftel signatures
kde_cons(GSEs, 'Nef')
# Verhaak signatures
kde_cons(GSEs, 'Ver')
# Mix 
kde_cons(GSEs, 'Mix')