#%%
import os
import glob
import itertools as it
import pandas as pd
import numpy    as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''
sns.set_context('paper', font_scale=1.8)
#%%
def corr_scraper(GSEs):
    df_mega = pd.DataFrame()
    # Iterate through GSEs
    for GSE in GSEs:
        # List GSMs in each GSE
        files = glob.glob('../Output/'+GSE+'/GSEA/*-corr.tsv')
        GSMs = [os.path.basename(file).replace('-corr.tsv','') for file in files]
        # Iterate through GSMs
        for GSM in GSMs:
            path = '../Output/'+GSE+'/GSEA/'
            # Read the correlation and p-value dataframes
            corr_df = pd.read_csv(path+GSM+'-corr.tsv', sep='\t', index_col=0)
            p_values = pd.read_csv(path+GSM+'-pval.tsv', sep='\t', index_col=0)
            # Convert matrix to long format
            corr_df = corr_df.stack().reset_index()
            p_values = p_values.stack().reset_index()
            # Merge the two dataframes
            df = pd.merge(corr_df, p_values, on=['level_0','level_1'])
            df.columns = ['src','tgt','corr','pval']
            # Add GSM and GSE columns
            df['GSM'] = GSM
            df['GSE'] = GSE
            # Append to mega dataframe
            df_mega = pd.concat([df_mega, df])
    return df_mega
#%%
def hmap_meta(df, sigs, tgt, pfx):
    # Select only the srcs with Signatures and the target
    df = df[df['src'].isin(sigs) & (df['tgt']==tgt)]
    # Pivot the dataframe
    df = df.pivot(index='GSM', columns='src', values='corr')
    # Reorder the columns
    df = df[sigs]
    # Plot the heatmap
    ax = sns.clustermap(data=df.T, cmap='coolwarm', vmax=1, vmin=-1, xticklabels=False, figsize=(8,6))
    ax.ax_heatmap.set_xlabel('')
    ax.ax_heatmap.set_ylabel('Signature')
    # Make the cbar smaller in left corner
    pos = ax.ax_cbar.get_position()
    ax.ax_cbar.set_position([pos.x0, pos.y0+0.05, pos.width, pos.height-0.05])
    plt.savefig('../figures/GSEA/Hmap_'+pfx+'_'+tgt+'.svg')
    # plt.close()
    # plt.clf()
#%%
def hmap_combi(df, sigs, title):
    # Make combination of src and tgt
    df['Combination'] = df['src'] + '-' + df['tgt']
    # Make combinations of signatures
    combi = list(it.combinations(sigs,2))
    combi = ['-'.join(x) for x in combi]
    # Select based on the pairs
    df = df[df['Combination'].isin(combi)]
     # Pivot the dataframe
    df = df.pivot(index='GSM', columns='Combination', values='corr')
    # Reorder the columns
    df = df[combi]
    # Plot the heatmap
    ax = sns.clustermap(data=df, cmap='coolwarm', vmax=1, vmin=-1, yticklabels=False)
    ax.ax_heatmap.set_ylabel('')
    plt.savefig('../figures/GSEA/Hmap_'+title+'.svg')
#%%
os.makedirs('../figures/GSEA/', exist_ok=True)
#%%
datasets = pd.read_csv('../Input/Datasets_Bulk.csv')
GSEs = datasets.GSE
#%%
Sigs = ['NefMES','NefAC','NefOPC','NefNPC','VerMES','VerCL','VerNL','VerPN']
tgtSig = ['PD-L1','GLYC','KEGG_CC']
#%%
df = corr_scraper(GSEs)
#%%
for tgt in tgtSig:
    hmap_meta(df, Sigs, tgt, 'Bulk')
# %%
hmap_combi(df, Sigs, 'Bulk')
# %%
datasets = pd.read_csv('../Input/Datasets_SC.csv')
GSEs = datasets.GSE
#%%
df = corr_scraper(GSEs)
#%%
for tgt in tgtSig:
    hmap_meta(df, Sigs, tgt, 'SC')
# %%
hmap_combi(df, Sigs, 'SC')