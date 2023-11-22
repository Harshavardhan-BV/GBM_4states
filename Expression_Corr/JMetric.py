import os
import glob
import argparse
import numpy as np
import pandas as pd
import itertools as it

# Read the GSE ID from command line
parser = argparse.ArgumentParser(
                    prog='J-Metric.py',
                    description='Calculates the J-metric of the correlation between the expression of genes in signatures')
parser.add_argument('GSE', type=str, nargs=1, help='GSE ID')
GSE = parser.parse_args().GSE[0]

def corr_consistency(GSE, GSM, sigs, suff):
    # Read the correlation matrix
    corr_df = pd.read_csv('../Output/'+GSE+'/Correlation/'+GSM+'_correlation.tsv', index_col=0, sep='\t')
    # Pairwise combinations of signatures
    combinats = list(it.combinations(sigs,2))
    cons = np.empty(len(combinats))
    for i in range(len(combinats)):
        combi = combinats[i]
        genes = np.empty(2,dtype='object')
        self_sum = np.empty(2)
        for j in range(2):
            gs = gbmgenes[combi[j]].dropna()
            # select only the genes that are in the correlation matrix
            gs = gs[gs.isin(corr_df.index) & gs.isin(corr_df.columns)]
            genes[j] = gs
            # subset the within group and across group correlations
            self = corr_df.loc[gs,gs]
            # calculate the mean within group correlation
            self_sum[j] = (np.tril(self).sum() - len(self)/2)/self.size
        cross = corr_df.loc[genes[0],genes[1]]
        # calculate the mean across group correlation
        cross_sum = cross.values.sum()/cross.size
        # Score the consistency
        cons[i] = (np.sum(self_sum) - cross_sum)/2
    df = pd.DataFrame(combinats,columns=['G1','G2'])
    df['Consistency'] = cons
    df.to_csv('../Output/'+GSE+'/Correlation/'+GSM+'_consistency_'+suff+'.tsv',sep='\t',index=False)

# List of GSM IDs
files = glob.glob('../Output/'+GSE+'/Correlation/*_correlation.tsv')
GSMs = [os.path.basename(x).replace(f'_correlation.tsv','') for x in files]
# Read the genes of signatures
gbmgenes = pd.read_csv("../Signatures/GBM_signatures.csv")

# Iterate over GSM samples
for GSM in GSMs:
    print(GSM)
    # Neftel signatures
    sigs = ['NefMES','NefAC','NefOPC','NefNPC']
    corr_consistency(GSE, GSM, sigs, 'Nef')
    # Verhaak signatures
    sigs = ['VerMES','VerCL','VerNL','VerPN']
    corr_consistency(GSE, GSM, sigs, 'Ver')
    # 2D signatures
    sigs = ['NefNPC','NefMES', 'VerPN', 'VerMES']
    corr_consistency(GSE, GSM, sigs, 'Mix')