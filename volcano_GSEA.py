import os
import glob
import itertools as it
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''

GSEs = ['CCLE','TCGA','GSE131928', 'GSE168004','GSE182109']

def volcplot(ax, df):
    df = df.copy()
    df['-log10(p-value)'] = -np.log10(df['pval']+1e-300)
    # Colour if pval<0.05 and red if corr>0.5, blue if corr<-0.5
    df['Correlation'] = np.where(df['pval']<0.05, np.where(df['corr']>0.3, '+ve', np.where(df['corr']<-0.3, '-ve', 'NA')), 'NA')
    sns.scatterplot(data=df, x='corr', y='-log10(p-value)', hue='Correlation', palette=['grey','tab:red','tab:blue'],hue_order=['NA','+ve','-ve'], ax=ax)
    # add lines for corr=0.5 and corr=-0.5 and pval=0.05
    ax.axvline(x=0.3, color='tab:red', linestyle='--')
    ax.axvline(x=-0.3, color='tab:blue', linestyle='--')
    ax.axhline(y=-np.log10(0.05), color='grey', linestyle='--')
    # remove legend
    ax.get_legend().remove()
    # # Add axis labels
    # ax.set_xlabel('Spearmans correlation')
    # ax.set_ylabel('-log10(p-value)')
    

def volc_grid(df, gs, suff):
    fig, ax = plt.subplots(2,3, figsize=(15,6), sharex=True, sharey=True)
    combi = list(it.combinations(gs, 2))
    for i in range(len(combi)):
        j = i//3
        k = i%3
        df_sub = df[(df['src']==combi[i][0]) & (df['tgt']==combi[i][1])]
        volcplot(ax[j,k], df_sub)
        ax[j,k].set_title(combi[i][0]+' vs '+combi[i][1])
    # Add common legend
    handles, labels = ax[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=3)
    plt.tight_layout()
    plt.savefig('./figures/GSEA/volcano_'+suff+'.svg')
    plt.clf()
    plt.close()

def corr_scraper(GSEs):
    df_mega = pd.DataFrame()
    # Iterate through GSEs
    for GSE in GSEs:
        # List GSMs in each GSE
        files = glob.glob('./Output/'+GSE+'/GSEA/*-corr.tsv')
        GSMs = [os.path.basename(file).replace('-corr.tsv','') for file in files]
        # Iterate through GSMs
        for GSM in GSMs:
            # Read the correlation and p-value dataframes
            corr_df = pd.read_csv('./Output/'+GSE+'/GSEA/'+GSM+'-corr.tsv', sep='\t', index_col=0)
            p_values = pd.read_csv('./Output/'+GSE+'/GSEA/'+GSM+'-pval.tsv', sep='\t', index_col=0)
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
            
os.makedirs('./figures/GSEA', exist_ok=True)
df_mega = corr_scraper(GSEs)

sigs = ['NefMES','NefAC','NefOPC','NefNPC']
volc_grid(df_mega, sigs, 'Nef')

sigs = ['VerMES','VerCL','VerNL','VerPN']
volc_grid(df_mega, sigs, 'Ver')

sigs = ['VerMES','NefMES','NefNPC','VerPN']
volc_grid(df_mega, sigs, '2D')