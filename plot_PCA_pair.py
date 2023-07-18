import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''

def loada_barplot(GSE, GSM, suff):
    print(GSM,suff)
    # Load the loadings matrix
    load_df = pd.read_csv('./Output/'+GSE+'/PCA/'+GSM+'_loadings_'+suff+'.tsv',sep='\t',index_col=0)
    # Load the variance explained
    var_df = pd.read_csv('./Output/'+GSE+'/PCA/'+GSM+'_expvar_'+suff+'.tsv',sep='\t',header=None)
    var_df = var_df*100
    # Read the GBM signature
    gbmgenes = pd.read_csv("./Signatures/GBM_signatures.csv")
    # Select only the columns that start with suff
    gs = suff.split('-')
    gbmgenes = gbmgenes.loc[:,gs]
    # select only the genes in the signature and merge with load_df
    gbmgenes = gbmgenes.melt(var_name='Signature').dropna()
    load_df = load_df.merge(gbmgenes,left_index=True,right_on='value', how='left')
    # Plot the loadings matrix but dont plot the ylabels
    g = sns.barplot(data=load_df,x='PC1',y=load_df.index, hue='Signature', dodge=False, orient='h', order= load_df.sort_values('PC1',ascending=False).index)
    g.set(yticklabels=[], yticks=[])
    plt.xlabel('PC1\n'+str(round(var_df.values[0][0],2))+'%')
    plt.tight_layout()
    plt.savefig('./figures/'+GSE+'/PCA/'+GSM+'_loadings_barplot_'+suff+'.svg')
    plt.close()
    plt.clf()

def expvar_plot(GSE,GSM,gs,suff):
    gs.append(suff)
    li = []
    # Read the variance explained
    for g in gs:
        df = pd.read_csv('./Output/'+GSE+'/PCA/'+GSM+'_expvar_'+g+'.tsv',sep='\t',header=None)
        df = df*100
        # Concatenate the variance explained
        li.append([df.iloc[0,0]])
    # Make a dataframe of the variance explained with index the gs
    df = pd.DataFrame(li,index=gs, columns=['Variance_explained'])
    sns.barplot(data=df, x=df.index, y='Variance_explained')
    plt.savefig('./figures/'+GSE+'/PCA/'+GSM+'_expvar_'+suff+'.svg')
    plt.close()
    plt.clf()
    
# GSE ID
# GSE = "GSE168004"
GSE = "GSE131928"
# GSE = "GSE182109"
# GSE = "CCLE"
# GSE = "TCGA"
# Make the output directory
os.makedirs('figures/'+GSE+'/PCA/', exist_ok=True)
# List files in Output
files = glob.glob('Output/'+GSE+'/PCA/*_loadings_Nef.tsv')
# Get the GSMs
GSMs = [os.path.basename(x).replace(f'_loadings_Nef.tsv','') for x in files]
# Iterate over GSM samples
for GSM in GSMs:
    gs = ['NefNPC-VerMES','NefOPC-NefAC']
    for suff in gs:
        loada_barplot(GSE,GSM,suff)
    expvar_plot(GSE,GSM, gs,'Nef')
    gs = ['VerPN-VerMES','VerCL-VerNL']
    for suff in gs:
        loada_barplot(GSE,GSM,suff)
    expvar_plot(GSE,GSM,gs,'Ver')