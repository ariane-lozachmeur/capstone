import pandas as pd
import numpy as np
from math import *
import matplotlib.pyplot as plt
import matplotlib
import matplotlib
import scipy.stats


def load_aggr(file,column):
    df_ces = pd.read_table(file)
    ces = df_ces.sum()[2:]
    return pd.DataFrame(ces,columns=[column])


def load_data(file, cat=None):    
    df = pd.read_table(file)
    patients = df['COMMON']
    df = df.transpose()[2:]
    df.columns = [s + '_' +cat for s in patients]
    return df

def load_mutation(file):
    df = pd.read_table(file)[['Sample ID', 'Mutation Count', 'CNA']]
    return df

from pylab import plot, show, savefig, xlim, figure, \
                hold, ylim, legend, boxplot, setp, axes

# function for setting the colors of the box plots pairs
def setBoxColors(bp):
    setp(bp['boxes'][0], color='blue')
    setp(bp['caps'][0], color='blue')
    setp(bp['caps'][1], color='blue')
    setp(bp['whiskers'][0], color='blue')
    setp(bp['whiskers'][1], color='blue')
    setp(bp['fliers'][0], color='blue')
    setp(bp['fliers'][1], color='blue')
    setp(bp['medians'][0], color='blue')

    setp(bp['boxes'][1], color='red')
    setp(bp['caps'][2], color='red')
    setp(bp['caps'][3], color='red')
    setp(bp['whiskers'][2], color='red')
    setp(bp['whiskers'][3], color='red')
    # setp(bp['fliers'][2], color='red')
    # setp(bp['fliers'][3], color='red')
    setp(bp['medians'][1], color='red')


def boxplotcomp(df, genes, selectedgene=None,):
    d = {}
    for g in genes:
        if selectedgene == None:
            d[g] = [list(df[df['msi_status']<0.74][g+'_expr']),list(df[df['msi_status']>0.74][g+'_expr'])]
        else:
            li = []
            for ix,row in df.iterrows():
                if selectedgene in row['msi_deficiency']:
                    li.append(True)
                else:
                    li.append(False)

            d[g] = [list(df[li][g+'_expr']),list(df[df['msi_status']==1][g+'_expr'])]
    figlen = 20
    figheigth = 12

    fig = figure(figsize=(figlen, figheigth))
    ax = axes()
    hold(True)

    i = 0
    for g in genes:
        bp = boxplot(d[g], positions=[i+1,i+2],widths=0.6)
        setBoxColors(bp)
        i+=2

    # set axes limits and labels
    plt.xticks(range(0, (len(genes)+1) * 2, 2), [''] + genes)
    if selectedgene == None:
        plt.title('MSS and MSI tumors comparison')

    else:
        plt.title(selectedgene+' deficient ('+str(len(d[selectedgene][0]))+' tumors) and MSI tumors comparison')
    # draw temporary red and blue lines and use them to create a legend
    hR, = plot([1,1],'r-')
    hB, = plot([1,1],'b-')
    if selectedgene == None:
        legend((hR, hB),('MSI tumors','MSS'))
    else:
        legend((hR, hB),('MSI tumors',selectedgene+' deficient tumors'))
    hB.set_visible(False)
    hR.set_visible(False)

    if selectedgene == None:
        savefig('results/MSS and MSI tumors comparison.png')
    else:
        savefig('results/'+selectedgene+' deficient and MSI tumors comparison.png')
    
    plt.show()


def get_mmr_status(df, threshs, genes_full, genes_part):
    df['mmr_status'] = np.zeros(len(df))
    df['mmr_deficiency'] = [list() for x in range(len(df.index))]
    for g in genes_full:
        # for each gene of this list, if the expression is under the threshold, we consider it to be dMMR.
        for ix in df[df[g+'_expr']< threshs[0] ]['mmr_status'].index:
            df = df.set_value(ix,'mmr_status',1)
            df.loc[ix,'mmr_deficiency'] == df.loc[ix,'mmr_deficiency'].append(g)

        # if the gene presents one mutation, we consider the tumor dMMR.
        no_mutations = pd.isnull(df[g+'_mut'])
        for ix in no_mutations[no_mutations==False].index:
            df = df.set_value(ix,'mmr_status',1)
            df.loc[ix,'mmr_deficiency'] == df.loc[ix,'mmr_deficiency'].append(g)

    for g in genes_part:      
        for ix in df[df[g+'_expr']<threshs[1] ]['mmr_status'].index:
            try:
                status = float(df[ix:ix+1]['mmr_status'])
            except:
                status = -1
            
            if status == 0.:
                df = df.set_value(ix,'mmr_status',0.25)
            elif status == 0.5:
                df = df.set_value(ix,'mmr_status',0.75)

    return df
