import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import itertools as it
plt.rcParams["svg.hashsalt"]=''

def corr_plot(GSE, GSM):
    # Read the correlation matrix
    corr_df = pd.read_csv('./Output/'+GSE+'/Correlation/'+GSM+'_correlation.tsv',sep='\t',index_col=0)
    # Remove the NaNs
    corr_df.dropna(axis=0, inplace=True)
    # Read the GBM genes
    gbmgenes = pd.read_excel("./Signatures/1-s2.0-S0092867419306877-mmc2.xlsx",skiprows=4)
    gbmgenes = gbmgenes.iloc[:,0:6]
    # make a colour map for the genes
    lut = {
        'MES1': 'tab:red',
        'MES2': 'tab:red',
        'NPC1': 'tab:blue',
        'NPC2': 'tab:blue',
        'OPC': 'tab:green',
        'AC': 'tab:orange',
    }
    # map the colours to the genes
    row_colors = corr_df.index.map(gbmgenes.melt().drop_duplicates(subset=['value']).set_index('value')['variable'])
    row_colors = row_colors.map(lut)
    row_colors = row_colors.fillna('white')
    # plot the clustermap of the correlation matrix
    sns.clustermap(data=corr_df, cmap='coolwarm',vmax=1, vmin=-1, row_colors=row_colors)
    # save the figure
    plt.savefig('./figures/'+GSE+'/Correlation/'+GSM+'_Corrplot_coloured.svg')
    plt.close()
    plt.clf()

def gene_subset(gbmgenes, state):
    if state == 'MES':
        return pd.concat([gbmgenes['MES1'].dropna(), gbmgenes['MES2'].dropna()],ignore_index=True)
    elif state == 'NPC':
        return pd.concat([gbmgenes['NPC1'].dropna(), gbmgenes['NPC2'].dropna()],ignore_index=True)
    else:
        return gbmgenes[state].dropna()

def corr_consistency(GSE, GSM):
    # Read the correlation matrix
    corr_df = pd.read_csv('./Output/'+GSE+'/Correlation/'+GSM+'_correlation.tsv',sep='\t',index_col=0)
    # Remove the NaNs
    corr_df.dropna(axis=0, inplace=True)
    # Read the GBM genes
    gbmgenes = pd.read_excel("./Signatures/1-s2.0-S0092867419306877-mmc2.xlsx",skiprows=4)
    gbmgenes = gbmgenes.iloc[:,0:6]
    # states to compare
    states = ['MES','AC','OPC','NPC']
    cons = []
    combinats = list(it.combinations(states,2))
    for combi in combinats:
        genes1 = gene_subset(gbmgenes, combi[0])
        genes2 = gene_subset(gbmgenes, combi[1])
        # select only the genes that are in the correlation matrix
        genes1 = genes1[genes1.isin(corr_df.index) & genes1.isin(corr_df.columns)]
        genes2 = genes2[genes2.isin(corr_df.index) & genes2.isin(corr_df.columns)]
        # subset the within group and across group correlations
        self1 = corr_df.loc[genes1,genes1]
        self2 = corr_df.loc[genes2,genes2]
        cross = corr_df.loc[genes1,genes2]
        # calculate the mean within group correlation
        self1_sum = (np.tril(self1).sum() - len(self1)/2)/self1.size
        self2_sum = (np.tril(self2).sum() - len(self2)/2)/self2.size
        # calculate the mean across group correlation
        cross_sum = cross.values.sum()/cross.size
        # Score the consistency
        consistency = (self1_sum + self2_sum - cross_sum)/2
        cons.append(consistency)
    df = pd.DataFrame(combinats,columns=['G1','G2'])
    df['Consistency'] = cons
    df.to_csv('./Output/'+GSE+'/Correlation/'+GSM+'_consistency.tsv',sep='\t',index=False)
    plt.bar(range(len(combinats)),cons)
    plt.xticks(range(len(combinats)), combinats, rotation=45)
    plt.ylabel('Consistency')
    plt.tight_layout()
    plt.savefig('./figures/'+GSE+'/Correlation/'+GSM+'_Consistency.svg')
    plt.clf()
    plt.close()
    

# GSE ID     
GSE = "GSE168004"
# List of GSM IDs
GSMs = [x.replace("_correlation.tsv",'') for x in os.listdir('Output/'+GSE+'/Correlation/')]
# Make directory for figures
os.makedirs('figures/'+GSE+'/Correlation/', exist_ok=True)

# Iterate over GSM samples
for GSM in GSMs:
    corr_plot(GSE, GSM)
    corr_consistency(GSE,GSM)
