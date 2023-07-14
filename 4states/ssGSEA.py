import os
import glob
import gseapy as gp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from multiprocessing import cpu_count

def plot_seprate(df, suff, GSE, GSM):
    sns.clustermap(df.corr(method='spearman'), cmap='coolwarm',vmax=1, vmin=-1)
    plt.savefig('figures/'+GSE+'/ssGSEA/ssGSEA_'+GSM+'_corrplot_'+suff+'.svg')
    plt.clf()
    plt.close()

def run_ssgsea(GSE, GSM):
    # Read the counts
    df = pd.read_csv("./Data_generated/"+GSE+"/Counts/"+GSM+"_counts.tsv",sep='\t', index_col=0)
    # Get the gene sets
    geneSets = './Signatures/GBM_states.gmt'
    # Calculate ssGSEA scores
    ss = gp.ssgsea(data=df,
               gene_sets=geneSets,
               outdir=None,
               sample_norm_method='rank', # choose 'custom' will only use the raw value of `data`
               no_plot=True, threads=cpu_count()-2)
    # Convert to dataframe with each signature as a column
    df = ss.res2d.pivot(index='Name', columns='Term', values='NES')
    df = df.convert_dtypes()
    # Save the output
    df.to_csv("./Output/"+GSE+"/ssGSEA/"+GSM+"-ssgsea.csv", index=True)
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
GSE = 'CCLE'
# GSE = 'TCGA'
# List files in Output
files = glob.glob('Data_generated/'+GSE+'/Counts/*_counts.tsv')
# Get the GSMs
GSMs = [os.path.basename(x).replace(f'_counts.tsv','') for x in files]
# Make directory for figures and output
os.makedirs('Output/'+GSE+'/ssGSEA/', exist_ok=True)
os.makedirs('figures/'+GSE+'/ssGSEA/', exist_ok=True)

for i in range(len(GSMs)):
    run_ssgsea(GSE, GSMs[i])