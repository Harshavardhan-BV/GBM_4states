import os
import glob
import itertools as it
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''

def volcplot(ax, df):
    df = df.copy()
    df['-log10(p-value)'] = -np.log10(df['pval'])
    # replace infinite pvalue to the max 
    max_pval = df['-log10(p-value)'].replace(np.inf, np.nan).max()
    df['-log10(p-value)'] = np.where(df['pval']==0, max_pval, df['-log10(p-value)'])
    # Colour if pval<0.05 and red if corr>0.5, blue if corr<-0.5
    df['Correlation'] = np.where(df['pval']<0.05, np.where(df['corr']>0.3, '+ve', np.where(df['corr']<-0.3, '-ve', 'NA')), 'NA')
    # Plot the scatter
    sns.scatterplot(data=df, x='corr', y='-log10(p-value)', hue='Correlation', palette=['grey','tab:red','tab:blue'],hue_order=['NA','+ve','-ve'], ax=ax)
    # Add counts of +ve and -ve
    ax.text(0.9, 0.9, '+ve: '+str(sum(df['Correlation']=='+ve')), transform=ax.transAxes, ha='right')
    ax.text(0.1, 0.9, '-ve: '+str(sum(df['Correlation']=='-ve')), transform=ax.transAxes, ha='left')
    ax.text(0.5, 0.9, 'NA: '+str(sum(df['Correlation']=='NA')), transform=ax.transAxes, ha='center')
    # add lines for corr=0.3 and corr=-0.3 and pval=0.05
    ax.axvline(x=0.3, color='tab:red', linestyle='--')
    ax.axvline(x=-0.3, color='tab:blue', linestyle='--')
    ax.axhline(y=-np.log10(0.05), color='grey', linestyle='--')
    # remove legend
    ax.get_legend().remove()
    

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
    # Change xlabels
    for i in range(3):
        ax[1,i].set_xlabel('Correlation')
    # Add the legend beyond the plot bound and resize the plot
    fig.legend(handles, labels, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.2))
    fig.subplots_adjust(bottom=0.2)
    fig.tight_layout()
    plt.savefig('../figures/GSEA/volcano_'+suff+'.svg')
    plt.clf()
    plt.close()

def volc_grid_2D(df, gtgt, gsrc, suff):
    ls = len(gsrc)
    lt = len(gtgt)
    fig, ax = plt.subplots(ls, lt, figsize=(5*lt,3*ls), sharex=True, sharey=True, squeeze=False)
    for j in range(lt):
        for k in range(ls):
            src = gsrc[k]
            tgt = gtgt[j]
            df_sub = df[(df['src']==src) & (df['tgt']==tgt)]
            volcplot(ax[k,j], df_sub)
            ax[k,j].set_title(src+' vs '+tgt)
    # Add common legend
    handles, labels = ax[0,0].get_legend_handles_labels()
    # Change xlabels
    for i in range(lt):
        ax[-1,i].set_xlabel('Correlation')
    # Add the legend beyond the plot bound and resize the plot
    fig.legend(handles, labels, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.2))
    fig.subplots_adjust(bottom=0.2)
    fig.tight_layout()
    plt.savefig('../figures/GSEA/volcano_'+suff+'.svg')
    plt.clf()
    plt.close()

def corr_scraper(GSEs, expr=False):
    df_mega = pd.DataFrame()
    # Iterate through GSEs
    for GSE in GSEs:
        # List GSMs in each GSE
        files = glob.glob('../Output/'+GSE+'/GSEA/*-corr.tsv')
        GSMs = [os.path.basename(file).replace('-corr.tsv','') for file in files]
        # Iterate through GSMs
        for GSM in GSMs:
            if expr:
                path = '../Output/'+GSE+'/GSEA-Expr/'
            else:
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
            
os.makedirs('../figures/GSEA', exist_ok=True)

datasets = pd.read_csv('../Input/Datasets_Bulk.csv')
GSEs = datasets.GSE

df_mega = corr_scraper(GSEs)
df_mega_expr = corr_scraper(GSEs, expr=True)

sigs_nef = ['NefMES','NefAC','NefOPC','NefNPC']
sigs_ver = ['VerMES','VerCL','VerNL','VerPN']
sigs_2D = ['VerMES','NefMES','NefNPC','VerPN']
sig_met = ['OXPHOS','GLYC','FAO']
sig_cc = ['G2M','KEGG_CC', 'WP_CC','BIOCARTA_CC']
sig_imm = ['PD-L1']
gene_imm = ['CTLA4','CD274','LAG3','HAVCR2','CD47','LGALS9','CD276']

volc_grid(df_mega, sigs_nef, 'Bulk_Nef')
volc_grid_2D(df_mega, sigs_nef, sig_met, 'Bulk_Nef_Met')
volc_grid_2D(df_mega, sigs_nef, sig_cc, 'Bulk_Nef_CC')
volc_grid_2D(df_mega, sigs_nef, sig_imm, 'Bulk_Nef_Imm')
volc_grid_2D(df_mega_expr, sigs_nef, gene_imm, 'Bulk_Nef_Imm_Expr')
volc_grid(df_mega, sigs_ver, 'Bulk_Ver')
volc_grid_2D(df_mega, sigs_ver, sig_met, 'Bulk_Ver_Met')
volc_grid_2D(df_mega, sigs_ver, sig_cc, 'Bulk_Ver_CC')
volc_grid_2D(df_mega, sigs_ver, sig_imm, 'Bulk_Ver_Imm')
volc_grid_2D(df_mega_expr,  sigs_ver, gene_imm, 'Bulk_Ver_Imm_Expr')
volc_grid(df_mega, sigs_2D, 'Bulk_2D')


datasets = pd.read_csv('../Input/Datasets_SC.csv')
GSEs = datasets.GSE

df_mega = corr_scraper(GSEs)
df_mega_expr = corr_scraper(GSEs, expr=True)

volc_grid(df_mega, sigs_nef, 'SC_Nef')
volc_grid_2D(df_mega, sigs_nef, sig_met, 'SC_Nef_Met')
volc_grid_2D(df_mega, sigs_nef, sig_cc, 'SC_Nef_CC')
volc_grid_2D(df_mega, sigs_nef, sig_imm, 'SC_Nef_Imm')
volc_grid_2D(df_mega_expr, sigs_nef, gene_imm, 'SC_Nef_Imm_Expr')
volc_grid(df_mega, sigs_ver, 'SC_Ver')
volc_grid_2D(df_mega, sigs_ver, sig_met, 'SC_Ver_Met')
volc_grid_2D(df_mega, sigs_ver, sig_cc, 'SC_Ver_CC')
volc_grid_2D(df_mega, sigs_ver, sig_imm, 'SC_Ver_Imm')
volc_grid_2D(df_mega_expr,  sigs_ver, gene_imm, 'SC_Ver_Imm_Expr')
volc_grid(df_mega, sigs_2D, 'SC_2D')