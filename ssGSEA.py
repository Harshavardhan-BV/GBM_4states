import os
import glob
import gseapy as gp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from multiprocessing import cpu_count

def run_ssgsea(GSE, GSM):
    # Read the counts
    df = pd.read_csv("./Data_generated/"+GSE+"/Counts/"+GSM+"_counts.tsv",sep='\t', index_col=0)
    # Get the gene sets
    geneSets = './Signatures/GBM_6states.gmt'
    # Calculate ssGSEA scores
    ss = gp.ssgsea(data=df,
               gene_sets=geneSets,
               outdir=None,
               sample_norm_method='rank', # choose 'custom' will only use the raw value of `data`
               no_plot=True, threads=cpu_count()-2)
    df = ss.res2d.pivot(index='Name', columns='Term', values='NES')
    df = df.convert_dtypes()
    df.to_csv("./Output/"+GSE+"/ssGSEA/"+GSM+"-ssgsea.csv", index=True)
    df = df.loc[:, (df != 0).any(axis=0)]
    fourdf = df[['NPC','AC','OPC','MES']]
    sns.clustermap(df.corr(method='spearman'), cmap='coolwarm',vmax=1, vmin=-1)
    plt.savefig('figures/'+GSE+'/ssGSEA/ssGSEA_'+GSM+'_corrplot_all.svg')
    plt.clf()
    plt.close()
    sns.clustermap(fourdf.corr(method='spearman'), cmap='coolwarm',vmax=1, vmin=-1)
    plt.savefig('figures/'+GSE+'/ssGSEA/ssGSEA_'+GSM+'_corrplot_4D.svg')
    plt.clf()
    plt.close()

# GSE ID
GSE = 'CCLE'
# List files in Output
files = glob.glob('Data_generated/'+GSE+'/Counts/*_counts.tsv')
# Get the GSMs
GSMs = [os.path.basename(x).replace(f'_counts.tsv','') for x in files]
# Make directory for figures and output
os.makedirs('Output/'+GSE+'/ssGSEA/', exist_ok=True)
os.makedirs('figures/'+GSE+'/ssGSEA/', exist_ok=True)

for i in range(len(GSMs)):
    run_ssgsea(GSE, GSMs[i])