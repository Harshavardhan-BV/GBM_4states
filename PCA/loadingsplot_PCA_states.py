import os
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''
sns.set_context('paper',font_scale=1.8)

# Read the GSE ID from command line
parser = argparse.ArgumentParser(
                    prog='loadingsplot_PCA_states.py',
                    description='Plots a screeplot and barplot of the loadings of the PCA')
parser.add_argument('GSE', type=str, nargs=1, help='GSE ID')
GSE = parser.parse_args().GSE[0]

def screeplot(GSE,GSM,suff):
    print(GSM,suff)
    # Load the variance explained
    var_df = pd.read_csv('./Output/'+GSE+'/PCA/'+GSM+'_expvar_'+suff+'.tsv',sep='\t',header=None)
    var_df = var_df*100
    # Select only first 10
    var_df = var_df.iloc[:10,:]
    # Plot the variance explained
    plt.figure(figsize=(5,5))
    plt.plot(var_df.index+1,var_df.values)
    plt.xlabel('Principal Component')
    plt.ylabel('Variance Explained (%)')
    plt.tight_layout()
    plt.savefig('./figures/'+GSE+'/PCA/'+GSM+'_screeplot_'+suff+'.svg')
    plt.close()
    plt.clf()

# Sort the loadings of PCi and plot a barplot. Colour the bars according to the geneset
def loada_barplot(GSE, GSM, PCi, suff):
    print(GSM,suff)
    # Load the loadings matrix
    load_df = pd.read_csv('../Output/'+GSE+'/PCA/'+GSM+'_loadings_'+suff+'.tsv',sep='\t',index_col=0)
    # Load the variance explained
    var_df = pd.read_csv('../Output/'+GSE+'/PCA/'+GSM+'_expvar_'+suff+'.tsv',sep='\t',header=None)
    var_df = var_df*100
    # Read the GBM signature
    gbmgenes = pd.read_csv("../Signatures/GBM_signatures.csv")
    # Select only the columns that start with suff
    gbmgenes = gbmgenes.filter(like=suff)
    # select only the genes in the signature and merge with load_df
    gbmgenes = gbmgenes.melt(var_name='Signature').dropna()
    load_df = load_df.merge(gbmgenes,left_index=True,right_on='value', how='left')
    # Set colour map as per the geneset
    lut = { col: f'tab:{clr}' for col, clr in zip(['NefMES','NefAC', 'NefOPC', 'NefNPC', 'VerMES', 'VerCL', 'VerNL', 'VerPN'],['red','orange','green','blue','red','orange','green','blue'])}
    # Plot the loadings matrix but dont plot the ylabels
    plt.figure(figsize=(7,5))
    g = sns.barplot(data=load_df,x='PC'+str(PCi),y=load_df.index, hue='Signature', dodge=False, orient='h', order= load_df.sort_values('PC'+str(PCi),ascending=False).index, palette=lut)
    g.set(yticklabels=[], yticks=[])
    plt.xlabel('PC'+str(i)+'\n'+str(round(var_df.values[PCi-1][0],2))+'%')
    plt.tight_layout()
    plt.savefig('../figures/'+GSE+'/PCA/'+GSM+'_loadings_PC'+str(PCi)+'_barplot_'+suff+'.svg')
    plt.close()
    plt.clf()

# Get the top 10% genes and plot a piechart of the proportions of genes in each signature
def topload(GSE,GSM,PCi,suff):
    print(GSM,suff)
    # Load the loadings matrix
    load_df = pd.read_csv('../Output/'+GSE+'/PCA/'+GSM+'_loadings_'+suff+'.tsv',sep='\t',index_col=0)
    # Read the GBM signature
    gbmgenes = pd.read_csv("../Signatures/GBM_signatures.csv")
    # Select only the columns that start with suff
    gbmgenes = gbmgenes.filter(like=suff)
    # select only the genes in the signature and merge with load_df
    gbmgenes = gbmgenes.melt(var_name='Signature').dropna()
    load_df = load_df.merge(gbmgenes,left_index=True,right_on='value', how='left')
    # Set colour map as per the geneset
    lut = { col: f'tab:{clr}' for col, clr in zip(['NefMES','NefAC', 'NefOPC', 'NefNPC', 'VerMES', 'VerCL', 'VerNL', 'VerPN'],['red','orange','green','blue','red','orange','green','blue'])}
    # Select only the top 10% genes +ve and bottom 10% -ve genes
    n = int(len(load_df)*0.1)
    topload_df = load_df.sort_values('PC'+str(PCi),ascending=False).iloc[:n,:]
    topload_df = topload_df[topload_df['PC'+str(PCi)]>0].groupby('Signature').size()
    botload_df = load_df.sort_values('PC'+str(PCi),ascending=True).iloc[:n,:]
    botload_df = botload_df[botload_df['PC'+str(PCi)]<0].groupby('Signature').size()
    # Check by absolute values as well
    load_df['PC'+str(PCi)] = load_df['PC'+str(PCi)].abs()
    load_df = load_df.sort_values('PC'+str(PCi),ascending=False).iloc[:n,:].groupby('Signature').size()
    # Plot a pie chart of proportions of genes in each signature
    plt.figure(figsize=(5,5))
    plt.pie(topload_df, labels=topload_df.index, colors=[lut[x] for x in topload_df.index])
    plt.savefig('../figures/'+GSE+'/PCA/'+GSM+'_PC'+str(PCi)+'_toppie_'+suff+'.svg')
    plt.close()
    plt.clf()
    plt.figure(figsize=(5,5))
    plt.pie(botload_df, labels=botload_df.index, colors=[lut[x] for x in botload_df.index])
    plt.savefig('../figures/'+GSE+'/PCA/'+GSM+'_PC'+str(PCi)+'_botpie_'+suff+'.svg')
    plt.close()
    plt.clf()

# Make the output directory
os.makedirs('../figures/'+GSE+'/PCA/', exist_ok=True)
# Iterative over the signatures
for suff in ['Nef','Ver']:
    # List files in Output
    files = glob.glob('../Output/'+GSE+'/PCA/*_loadings_'+suff+'.tsv')
    # Get the GSMs
    GSMs = [os.path.basename(x).replace(f'_loadings_'+suff+'.tsv','') for x in files]
    # Iterate over GSM samples
    for GSM in GSMs:
        screeplot(GSE,GSM, suff)
        for i in range(1,3):
            loada_barplot(GSE,GSM,i,suff)
            topload(GSE,GSM,i,suff)
