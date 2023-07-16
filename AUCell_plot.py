import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
plt.rcParams["svg.hashsalt"]=''

def spearman_pval(x,y):
    return spearmanr(x,y)[1]

def plot_seprate(df, suff, GSE, GSM):
    # Find the correlation of the matrix
    corr_df = df.corr(method='spearman')
    # Get the pvalues
    p_values = df.corr(method=spearman_pval)
    # change diagonal to 0. df.corr defaults diagonal to 1 for some reason
    np.fill_diagonal(p_values.values, 0)
    # Convert it to * format
    p_values = p_values.applymap(lambda x: '*' if x < 0.05 else '')
    # Plot the clustermap with the text being pvalues
    sns.clustermap(corr_df, cmap='coolwarm', vmax=1, vmin=-1, annot=p_values, fmt='', annot_kws={"size": 20})
    plt.savefig('figures/'+GSE+'/AUCell/AUCell_'+GSM+'_corrplot_'+suff+'.svg')
    plt.clf()
    plt.close()

def corr_plot_AUCell(GSE,GSM):
    print(GSM)
    # Read the AUCell scores
    df = pd.read_csv('Output/'+GSE+'/AUCell/'+GSM+'-AUCell.csv', index_col=0)
    df = df.loc[:, (df != 0).any(axis=0)]
    # Plot all signatures
    plot_seprate(df, 'all', GSE, GSM)
    # Select only Verhaak signatures
    df1 = df.loc[:, ['VerCL','VerMES','VerPN', 'VerNL']]
    plot_seprate(df1, 'Ver', GSE, GSM)
    # Select only Neftel signatures
    df1 = df.loc[:, ['NefMES','NefOPC','NefAC', 'NefNPC']]
    plot_seprate(df1, 'Nef', GSE, GSM)
    # Select only MES-MES vs NPC-PN
    df1 = df.loc[:, ['NefMES','VerMES','VerPN', 'NefNPC']]
    plot_seprate(df1, '2D', GSE, GSM)
    
     
# GSE ID
# GSE = "GSE168004"
# GSE = "GSE131928"
GSE = "GSE182109"

# List of GSM IDs
GSMs = [x.replace("-AUCell.csv",'') for x in os.listdir('Output/'+GSE+'/AUCell/')]
# Make directory for figures
os.makedirs('figures/'+GSE+'/AUCell/', exist_ok=True)

# Iterate over GSM samples
for GSM in GSMs:
    corr_plot_AUCell(GSE,GSM)
