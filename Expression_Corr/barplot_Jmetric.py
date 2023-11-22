import os
import glob
import argparse
import itertools as it
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''
sns.set_context('poster', font_scale=1.25)

# Read the GSE ID from command line
parser = argparse.ArgumentParser(
                    prog='barplot_Jmetric.py',
                    description='Plots barplot of the J-metric for each signature set')
parser.add_argument('GSE', type=str, nargs=1, help='GSE ID')
GSE = parser.parse_args().GSE[0]

def plot_consistency(GSE, GSM, suff):
    df = pd.read_csv('../Output/'+GSE+'/Correlation/'+GSM+'_consistency_'+suff+'.tsv',sep='\t')
    # Combine the columns G1 and G2
    df['Combination'] = df['G1'] + '-' + df['G2']
    plt.figure(figsize=(8,7))
    sns.barplot(data=df, x='Combination', y='Consistency')
    plt.ylabel('J-Metric')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig('../figures/'+GSE+'/Correlation/'+GSM+'_Consistency_'+suff+'.svg')
    plt.clf()
    plt.close()

# List of GSM IDs
files = glob.glob('../Output/'+GSE+'/Correlation/*_correlation.tsv')
GSMs = [os.path.basename(x).replace(f'_correlation.tsv','') for x in files]
# Make directory for figures
os.makedirs('../figures/'+GSE+'/Correlation/', exist_ok=True)

# Iterate over GSM samples
for GSM in GSMs:
    plot_consistency(GSE,GSM,'Nef')
    plot_consistency(GSE,GSM,'Ver')