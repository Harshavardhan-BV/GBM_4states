import os
import glob
import itertools as it
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''
sns.set_context('paper', font_scale=1.8)

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

def posneg_corr(df):
    df = df.copy()
    # Make correlation as 1 if + ve significant, -1 if -ve significant, 0 if not significant
    df['Correlation'] = np.where(df['pval']<0.05, np.where(df['corr']>0.3, '+ve', np.where(df['corr']<-0.3, '-ve', 'NA')), 'NA')
    # Group by src and tgt and count the number of +ve, -ve and 0 as separate columns
    df = df.groupby(['src','tgt','Correlation']).size().reset_index()
    df.columns = ['src','tgt','Correlation','Count']
    df = df.pivot(index=['src','tgt'], columns='Correlation', values='Count').reset_index()
    return df

def stacked_bar(df,sigs,tgt, pfx, signam):
    # Select only the srcs with Signatures and the target
    df = df[df['src'].isin(sigs) & (df['tgt']==tgt)]
    # Order as per sigs
    ys = ['NA','-ve','+ve']
    colors = ['grey','tab:blue','tab:red']
    # Plot the stacked barplot
    df.plot.bar(x='src',y=ys,color=colors,stacked=True, figsize=(6,6))
    plt.title(signam+' vs '+tgt)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.ylabel('Datasets')
    plt.xlabel('Signature')
    plt.tight_layout()
    plt.legend(loc='lower right')
    # Save the figure
    plt.savefig('../figures/GSEA/Bar_'+pfx+'_'+signam+'-'+tgt+'.svg')
    plt.close()
    plt.clf()

def stacked_bar_combi(df, sigs, title):
    df = df.copy()
    # Make combination of src and tgt
    df['Combination'] = df['src'] + '-' + df['tgt']
    # Make combinations of signatures
    combi = list(it.combinations(sigs,2))
    combi = ['-'.join(x) for x in combi]
    # Select based on the pairs
    df = df[df['Combination'].isin(combi)]
    # Order as per sigs
    ys = ['NA','-ve','+ve']
    colors = ['grey','tab:blue','tab:red']
    # Plot the stacked barplot
    df.plot.bar(x='Combination',y=ys,color=colors,stacked=True)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.ylabel('Datasets')
    # plt.xlabel('Combination')
    plt.tight_layout()
    plt.legend(loc='lower right')
    # Save the figure
    plt.savefig('../figures/GSEA/Bar_'+title+'.svg')
    plt.close()
    plt.clf()

os.makedirs('../figures/GSEA', exist_ok=True)

# Signatures
NefSig = ['NefMES','NefAC','NefOPC','NefNPC']
VerSig = ['VerMES','VerCL','VerNL','VerPN']
tgtSig = ['PD-L1','GLYC','KEGG_CC']
datasets = pd.read_csv('../Input/Datasets_Bulk.csv')
GSEs = datasets.GSE

df_count = corr_scraper(GSEs)
df_count = posneg_corr(df_count)

# Plot the stacked barplots for combinations
stacked_bar_combi(df_count, NefSig, 'Bulk_Nef')
stacked_bar_combi(df_count, VerSig, 'Bulk_Ver')
# Plot the stacked barplots
for tgt in tgtSig:
    stacked_bar(df_count,NefSig,tgt, 'Bulk', 'Nef')
    stacked_bar(df_count,VerSig,tgt, 'Bulk', 'Ver')

datasets = pd.read_csv('../Input/Datasets_SC.csv')
GSEs = datasets.GSE

df_count = corr_scraper(GSEs)
df_count = posneg_corr(df_count)

# Plot the stacked barplots for combinations
stacked_bar_combi(df_count, NefSig, 'SC_Nef')
stacked_bar_combi(df_count, VerSig, 'SC_Ver')
# Plot the stacked barplots
for tgt in tgtSig:
    stacked_bar(df_count,NefSig,tgt, 'SC', 'Nef')
    stacked_bar(df_count,VerSig,tgt, 'SC', 'Ver')