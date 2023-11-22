import os
import glob
import argparse
import pandas as pd
import seaborn as sns
import itertools as it
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

def aucorr(GSE, GSM, suff):
    print(GSM)
    df = pd.read_csv('../Output/'+GSE+'/'+score+'/'+GSM+'-'+score+'.csv', index_col=0)
    # Select columns with Nef in name
    df = df.filter(regex=suff)
    # Take combinations of columns
    combs = list(it.combinations(df.columns,2))
    # Assign cell type based on which score is max
    df['CellType'] = df.idxmax(axis=1)
    # Create colour paletter
    palette = {'NefNPC':'tab:blue', 'NefOPC':'tab:green', 'NefAC':'tab:orange', 'NefMES':'tab:red', 'VerPN':'tab:blue', 'VerNL':'tab:green', 'VerCL':'tab:orange', 'VerMES':'tab:red'}
    for combi in combs:
        # Plot the scatter plot
        sns.scatterplot(data=df, x=combi[0], y=combi[1], hue='CellType', edgecolor='none', s=10, palette=palette )
        plt.tight_layout()
        plt.savefig('../figures/'+GSE+'/GSEA/'+GSM+'-2D-{0}vs{1}.svg'.format(*combi))
        plt.close()
        plt.clf()

# List of GSM IDs
score = 'AUCell' if sc else 'ssGSEA'
files = glob.glob('../Output/'+GSE+'/'+score+'/*-'+score+'.csv')
GSMs = [os.path.basename(x).replace(f'-'+score+'.csv','') for x in files]
# Make directory for output
os.makedirs('../Output/'+GSE+'/GSEA/', exist_ok=True)

# Iterate over GSM samples
for GSM in GSMs:
    aucorr(GSE,GSM,'Nef')
    aucorr(GSE,GSM,'Ver')