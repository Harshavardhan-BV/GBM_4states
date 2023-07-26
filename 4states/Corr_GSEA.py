import os
import glob
import argparse
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# Read the GSE ID and sc from command line
parser = argparse.ArgumentParser(
                    prog='Corr_GSEA.py',
                    description='Calculates correlation and pvalues between the AUCell/ssGSEA scores for a given GSE')
parser.add_argument('GSE', type=str, nargs=1, help='GSE ID')
parser.add_argument('-sc', action=argparse.BooleanOptionalAction, help='Add this flag if single cell data')
GSE = parser.parse_args().GSE[0]
sc = parser.parse_args().sc

def spearman_pval(x,y):
    return spearmanr(x,y)[1]

def aucorr(GSE, GSM):
    print(GSM)
    df = pd.read_csv('Output/'+GSE+'/'+suff+'/'+GSM+'-'+suff+'.csv', index_col=0)
    # Find the correlation of the matrix
    corr_df = df.corr(method='spearman')
    # Get the pvalues
    p_values = df.corr(method=spearman_pval)
    # change diagonal to 0. df.corr defaults diagonal to 1 for some reason
    np.fill_diagonal(p_values.values, 0)
    # Save the correlation and pvalue matrices
    corr_df.to_csv('Output/'+GSE+'/GSEA/'+GSM+'-corr.tsv', sep='\t')
    p_values.to_csv('Output/'+GSE+'/GSEA/'+GSM+'-pval.tsv', sep='\t')

# List of GSM IDs
suff = 'AUCell' if sc else 'ssGSEA'
files = glob.glob('Output/'+GSE+'/'+suff+'/*-'+suff+'.csv')
GSMs = [os.path.basename(x).replace(f'-'+suff+'.csv','') for x in files]
# Make directory for output
os.makedirs('Output/'+GSE+'/GSEA/', exist_ok=True)

# Iterate over GSM samples
for GSM in GSMs:
    aucorr(GSE,GSM)