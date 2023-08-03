import os
import glob
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams["svg.hashsalt"]=''

# Read the GSE ID from command line
parser = argparse.ArgumentParser(
                    prog='2D_GSEA.py',
                    description='Plots a scatter plot of the PN/NPC vs MES scores')
parser.add_argument('GSE', type=str, nargs=1, help='GSE ID')
parser.add_argument('-sc', action=argparse.BooleanOptionalAction, help='Add this flag if single cell data')
GSE = parser.parse_args().GSE[0]
sc = parser.parse_args().sc

def aucorr(GSE, GSM):
    print(GSM)
    df = pd.read_csv('Output/'+GSE+'/'+suff+'/'+GSM+'-'+suff+'.csv', index_col=0)
    # Assign cell type based on which score is max
    df['CellType'] = df[['NefNPC','NefOPC','NefAC','NefMES']].idxmax(axis=1)
    # Create colour paletter
    palette = {'NefNPC':'tab:blue', 'NefOPC':'tab:green', 'NefAC':'tab:orange', 'NefMES':'tab:red', 'VerPN':'tab:blue', 'VerNL':'tab:green', 'VerCL':'tab:orange', 'VerMES':'tab:red'}
    # Plot the scatter plot
    sns.scatterplot(data=df, x='NefNPC', y='NefMES', hue='CellType', edgecolor='none', s=10, palette=palette )
    plt.tight_layout()
    plt.savefig('figures/'+GSE+'/GSEA/'+GSM+'-2D-Nef.svg')
    plt.close()
    plt.clf()
    # Assign cell type based on which score is max
    df['CellType'] = df[['VerPN','VerNL', 'VerCL', 'VerMES']].idxmax(axis=1)
    # Plot the scatter plot
    sns.scatterplot(data=df, x='VerPN', y='VerMES', hue='CellType',edgecolor='none', s=10, palette=palette)
    plt.tight_layout()
    plt.savefig('figures/'+GSE+'/GSEA/'+GSM+'-2D-Ver.svg')
    plt.close()
    plt.clf()

# List of GSM IDs
suff = 'AUCell' if sc else 'ssGSEA'
files = glob.glob('Output/'+GSE+'/'+suff+'/*-'+suff+'.csv')
GSMs = [os.path.basename(x).replace(f'-'+suff+'.csv','') for x in files]
# Make directory for output
os.makedirs('Output/'+GSE+'/GSEA/', exist_ok=True)

# Iterate over GSM samples
for GSM in GSMs:
    aucorr(GSE,GSM)