import os
import glob
import argparse
import numpy as np
import gseapy as gp
import pandas as pd
from multiprocessing import cpu_count

# Read the GSE ID from command line
parser = argparse.ArgumentParser(
                    prog='ssGSEA.py',
                    description='Calculates ssGSEA scores for a given GSE')
parser.add_argument('GSE', type=str, nargs=1, help='GSE ID')
GSE = parser.parse_args().GSE[0]

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
    df.to_csv("./Output/"+GSE+"/ssGSEA/"+GSM+"-ssGSEA.csv", index=True)

# List files in Output
files = glob.glob('Data_generated/'+GSE+'/Counts/*_counts.tsv')
# Get the GSMs
GSMs = [os.path.basename(x).replace(f'_counts.tsv','') for x in files]
# Make directory for figures and output
os.makedirs('Output/'+GSE+'/ssGSEA/', exist_ok=True)

for i in range(len(GSMs)):
    run_ssgsea(GSE, GSMs[i])