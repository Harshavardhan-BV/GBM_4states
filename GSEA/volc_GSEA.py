import os
import glob
import itertools as it
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''
sns.set_context('paper', font_scale=1.8)

def volcplot(df, src, tgt, pfx):
    df = df[(df['src']==src) & (df['tgt']==tgt)].copy()
    df['-log10(p-value)'] = -np.log10(df['pval'])
    # replace infinite pvalue to the max 
    max_pval = df['-log10(p-value)'].replace(np.inf, np.nan).max()
    df['-log10(p-value)'] = np.where(df['pval']==0, max_pval, df['-log10(p-value)'])
    # Colour if pval<0.05 and red if corr>0.5, blue if corr<-0.5
    df['Correlation'] = np.where(df['pval']<0.05, np.where(df['corr']>0.3, '+ve', np.where(df['corr']<-0.3, '-ve', 'NA')), 'NA')
    # Plot the scatter
    fig,ax = plt.subplots(figsize=(6,4.5))
    sns.scatterplot(data=df, x='corr', y='-log10(p-value)', hue='Correlation', palette=['grey','tab:red','tab:blue'],hue_order=['NA','+ve','-ve'], ax=ax)
    # Add counts of +ve and -ve
    ax.text(0.9, 0.9, '+ve: '+str(sum(df['Correlation']=='+ve')), transform=ax.transAxes, ha='right')
    ax.text(0.1, 0.9, '-ve: '+str(sum(df['Correlation']=='-ve')), transform=ax.transAxes, ha='left')
    ax.text(0.5, 0.9, 'NA: '+str(sum(df['Correlation']=='NA')), transform=ax.transAxes, ha='center')
    # add lines for corr=0.3 and corr=-0.3 and pval=0.05
    ax.axvline(x=0.3, color='tab:red', linestyle='--')
    ax.axvline(x=-0.3, color='tab:blue', linestyle='--')
    ax.axhline(y=-np.log10(0.05), color='grey', linestyle='--')
    ax.set_xlim(-1,1)
    # remove legend
    ax.get_legend().remove()
    ax.set_xlabel('Correlation')
    plt.title(src+" vs "+tgt)
    plt.tight_layout()
    # Save the figure
    plt.savefig('../figures/GSEA/Volcano_'+pfx+'_'+src+'-'+tgt+'.svg')
    plt.close()
    plt.clf()
    
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


# Signatures
NefSig = ['NefMES','NefAC','NefOPC','NefNPC']
VerSig = ['VerMES','VerCL','VerNL','VerPN']
tgtSig = ['PD-L1','GLYC','KEGG_CC']
datasets = pd.read_csv('../Input/Datasets_Bulk.csv')
GSEs = datasets.GSE

df_mega = corr_scraper(GSEs)
# Plot the volcano plots
volcplot(df_mega, 'NefMES', 'NefNPC', 'Bulk')
volcplot(df_mega, 'VerMES', 'VerPN', 'Bulk')
# Plot the volcano plots for
for src in ['NefMES','NefNPC','VerMES','VerPN']:
    for tgt in tgtSig:
        volcplot(df_mega,src,tgt,'Bulk')

datasets = pd.read_csv('../Input/Datasets_SC.csv')
GSEs = datasets.GSE

df_mega = corr_scraper(GSEs)
# Plot the volcano plots
volcplot(df_mega, 'NefMES', 'NefNPC', pfx='SC')
volcplot(df_mega, 'VerMES', 'VerPN', pfx='SC')
# Plot the volcano plots for
for src in ['NefMES','NefNPC','VerMES','VerPN']:
    for tgt in tgtSig:
        volcplot(df_mega,src,tgt,'SC')

