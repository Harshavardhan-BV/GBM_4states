import os
import glob
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams["svg.hashsalt"]=''

# Read the GSE ID from command line
parser = argparse.ArgumentParser(
                    prog='hmap_GSEA.py',
                    description='Plots a heatmap of the correlation between AUCell/ssGSEA scores')
parser.add_argument('GSE', type=str, nargs=1, help='GSE ID')
GSE = parser.parse_args().GSE[0]

def plot_seprate(corr_df,p_values, suff, GSE, GSM, sigs=None):
    if sigs:
        corr_df = corr_df.loc[sigs,sigs]
        p_values = p_values.loc[sigs,sigs]
    # Convert pvalues to * format
    p_values = p_values.applymap(lambda x: '*' if x < 0.05 else '')
    # Plot the clustermap with the text being pvalues
    sns.clustermap(corr_df, cmap='coolwarm', vmax=1, vmin=-1, annot=p_values, fmt='', annot_kws={"size": 20})
    plt.savefig('figures/'+GSE+'/GSEA/GSEA_'+GSM+'_corrplot_'+suff+'.svg')
    plt.clf()
    plt.close()

def corr_plot_AUCell(GSE,GSM):
    print(GSM)
    # Read the correlation and pvalue matrices
    corr_df = pd.read_csv('Output/'+GSE+'/GSEA/'+GSM+'-corr.tsv', index_col=0, sep='\t')
    p_values = pd.read_csv('Output/'+GSE+'/GSEA/'+GSM+'-pval.tsv', index_col=0, sep='\t')
    # Plot all signatures
    plot_seprate(corr_df, p_values, 'all', GSE, GSM)
    # Select only Verhaak signatures
    sigs = ['VerCL','VerMES','VerPN', 'VerNL']
    plot_seprate(corr_df, p_values, 'Ver', GSE, GSM, sigs=sigs)
    # Select only Neftel signatures
    sigs = ['NefMES','NefOPC','NefAC', 'NefNPC']
    plot_seprate(corr_df, p_values, 'Nef', GSE, GSM, sigs=sigs)
    # Select only MES-MES vs NPC-PN
    sigs = ['NefMES','VerMES','VerPN', 'NefNPC']
    plot_seprate(corr_df, p_values, '2D', GSE, GSM, sigs=sigs)
    
# List of GSM IDs
files = glob.glob('Output/'+GSE+'/GSEA/*-corr.tsv')
GSMs = [os.path.basename(x).replace(f'-corr.tsv','') for x in files]
# Make directory for figures
os.makedirs('figures/'+GSE+'/GSEA/', exist_ok=True)

# Iterate over GSM samples
for GSM in GSMs:
    corr_plot_AUCell(GSE,GSM)
