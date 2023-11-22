
import os
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''
sns.set_context('paper', font_scale=1.5)

# Read the GSE ID from command line
parser = argparse.ArgumentParser(
                    prog='Ninjastar_GSEA.py',
                    description='Makes a plot similar to Neftel et al (aka Ninja Star)  with AUCell scores')
parser.add_argument('GSE', type=str, nargs=1, help='GSE ID')
GSE = parser.parse_args().GSE[0]

def auninja(GSE, GSM):
    print(GSM)
    df = pd.read_csv('../Output/'+GSE+'/'+suff+'/'+GSM+'-'+suff+'.csv', index_col=0)
    # Define y axis
    df['y'] = df.loc[:,['NefNPC','NefOPC']].max(axis=1) - df.loc[:,['NefMES','NefAC']].max(axis=1)
    # Define x axis
    df['x'] = np.where(df['y']>0, (df.loc[:,'NefNPC']-df.loc[:,'NefOPC']), (df.loc[:,'NefMES']-df.loc[:,'NefAC']))
    # Assign cell type by argmax
    df['CellType'] = np.argmax(df.loc[:,['NefNPC','NefOPC','NefAC','NefMES']].values, axis=1)
    # Plot the scatter
    sns.scatterplot(data=df, x='x', y='y', hue='CellType', palette=['tab:blue','tab:green','tab:orange','tab:red'], size=1, linewidth=0, alpha=0.5)
    # Remove xlabel and ylabel and legend
    plt.xlabel('Relative metamodule score')
    plt.ylabel('Relative metamodule score')
    plt.legend([],[], frameon=False)
    # Set each corners as NPC, OPC, AC and MES
    plt.text(0.05, 0.9, 'OPC', transform=plt.gca().transAxes, ha='left', color='tab:green')
    plt.text(0.05, 0.05, 'AC', transform=plt.gca().transAxes, ha='left', color='tab:orange')
    plt.text(0.95, 0.05, 'MES', transform=plt.gca().transAxes, ha='right', color='tab:red')
    plt.text(0.95, 0.9, 'NPC', transform=plt.gca().transAxes, ha='right', color='tab:blue')
    # Put 0,0 in the middle
    y_max = max(abs(plt.ylim()[0]), abs(plt.ylim()[1]))
    x_max = max(abs(plt.xlim()[0]), abs(plt.xlim()[1]))
    plt.xlim(-x_max,x_max)
    plt.ylim(-y_max,y_max)
    plt.tight_layout()
    # Save the figure
    plt.savefig('../figures/'+GSE+'/GSEA/'+GSM+'-ninjastar.svg')
    # Close the figure
    plt.close()
    plt.clf()

# List of GSM IDs
suff = 'AUCell' 
files = glob.glob('../Output/'+GSE+'/'+suff+'/*-'+suff+'.csv')
GSMs = [os.path.basename(x).replace(f'-'+suff+'.csv','') for x in files]
# Make directory for output
os.makedirs('../figures/'+GSE+'/GSEA/', exist_ok=True)

# Iterate over GSM samples
for GSM in GSMs:
    auninja(GSE,GSM)
