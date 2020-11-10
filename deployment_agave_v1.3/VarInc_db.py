
import re
import pandas as pd
import os
import numpy as np
import itertools as it
import collections
from scipy.stats import norm


def recursively_default_dict():
    return collections.defaultdict(recursively_default_dict)



import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--db_dir', type=str,
                    default='.')

parser.add_argument('--species', type=str,
                    default='chimp')

parser.add_argument('--data', type=str,
                    default='sims')

parser.add_argument('--minS', type=int, # ignore pops if Nsamp < minS
                    default=5)

parser.add_argument('--samp', type=int, # if stepup= increment = max Nsamp, else = min (Nsamp).
                    default=200)

parser.add_argument('--steps', type=int, # Nber of Nsamp steps
                    default=50)

parser.add_argument('--stepup', type=str, # type of analysis - increment / other.
                    default='increment')

parser.add_argument('--simsN', type=int, # Numberr of simulations to read from available. random. 
                    default=0)           # if 0 uses all sims in sims_dir. 

parser.add_argument('--collapsed', type=bool, # collapse mutation counts.
                    default=False)

parser.add_argument('--pval', type=float, # Inclusion into H0 pval. 
                    default=0.01)           


args = parser.parse_args()


###############
species= args.species
############### plots

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from INFO_db import INFO_dict

########
########
fig_dir= 'Figures/VarSamp'
os.makedirs(fig_dir, exist_ok=True)
fig_dir= fig_dir + '/'


## directories
main_dir= os.getcwd() + '/'

pop_names_dict= {g:z for z,g in INFO_dict[species]['pop_dict'].items()}

######################################
######################################
##### READ DATA 
#######################################
from tools.dbcount import (
	info_array_collect
    )


db_dir= INFO_dict[species]['dirs'][args.data]
db_neigh= os.listdir(db_dir)


info_array, muts, counts= info_array_collect(db_dir)

#print(info_array[:10])

pops= INFO_dict[args.species]['pop_dict'].values()
pop_list= list(pops)

#####################################
########## III. COUNT COMPARISONS
##########
from tools.dbcount_cofactors import (
    sim_VarSub
)

stats_dict, counts_dict, ref_dict= sim_VarSub(counts,info_array,
                        row= row,col= col, si_max= args.samp)

#####################################
####################################
######## IV. Extract and Sort.
######## std distribution.
######## Across sub_samples:

from scipy.stats import norm
pval_threshold=  args.pval
std_threshold= norm.ppf(1-pval_threshold/2,loc= 0,scale= 1)

stat_used= ['self','ref']


pop_thresh_dict= {pop: {} for pop in list(pop_names_dict.values())}

for stat_request in stat_used:

    pop_diff_stats= {
        z: [stats_dict[z][x][stat_request] for x in sorted([x for x in stats_dict[z].keys() if x not in ['sizes', 'counts']])] for z in stats_dict.keys()
    }

    pop_diff_dict= {
        z: {
            'mean': [np.nanmean(x) for x in g],
            'std': [np.nanstd(x) for x in g],
        } for z,g in pop_diff_stats.items()
    }

    ######### Across references.
    ref_stats= {
        z: {
            'mean': np.nanmean(g),
            'ucl': np.nanmean(g) + std_threshold*np.nanstd(g),
            'lcl': np.nanmean(g) - std_threshold*np.nanstd(g),
            'std': np.nanstd(g)
        } for z,g in ref_dict.items()
    }

    average_cdf= {
        z: [norm.cdf(x,loc= ref_stats[z]['mean'],scale= ref_stats[z]['std']) for x in pop_diff_dict[z]['mean']] for z in pop_diff_dict.keys()
    }

    ####################################
    ####################################
    ########## IV. PLOT.
    ##########
    ########## i. ref_stats distribution - STD acrss sample sizes and for full-pops. 
    pop_colors= INFO_dict[args.species]['pop_colors']
    pop_names_dict= {g:z for z,g in INFO_dict[args.species]['pop_dict'].items()}

    from scipy.stats import norm

    ylab= 'sum of squared differences'
    xlab= 'sample size'
    for pop in pop_diff_dict.keys():
        pop_name= pop_names_dict[pop]
        surface= counts_dict[pop]['sizes']
        for plot_var in [0,1]:
            plt.figure(figsize=(10,10))
            
            
            y= pop_diff_dict[pop]['mean']
            error= pop_diff_dict[pop]['std']
            if plot_var:
                plt.errorbar(surface,y,yerr=error,label= 'var')   
            else:
                plt.plot(surface,y,label= 'mean var')

            line_surf= [min(counts_dict[pop]['sizes']),max(counts_dict[pop]['sizes'])]

            for z in ['mean','ucl','lcl']:
                y= [ref_stats[pop][z]] * 2
                plt.plot(line_surf,y,label= z)

            plt.xlabel(xlab)

            plt.ylabel(ylab)
            plt.title('grid SSD / sample proportion - control - {}'.format(pop_name), fontsize= 20)

            plt.legend()
            plt.savefig(fig_dir + sims_target + '_{}_gridSSD_{}_{}_sbsD.png'.format(stat_request,pop_name,["",'wSD'][plot_var]),bbox_inches='tight')
            plt.close()

    ### same plot.

    for set_ylim in [0,1]:
        plt.figure(figsize=(10,10))

        border_cols= {
            'mean': 'darkgrey',
            'ucl': 'darkgrey',
            'lcl': 'silver'
        }

        bars_dict= {}

        for pop in pop_diff_dict.keys():
            pop_name= pop_names_dict[pop]
            pop_col= pop_colors[pop_name]
            style= INFO_dict[args.species]['pop_style'][pop_name]
            surface= counts_dict[pop]['sizes']
            for plot_var in [0]:
                
                y= pop_diff_dict[pop]['mean']
                error= pop_diff_dict[pop]['std']

                vertical= [surface[x] for x in range(len(y)) if y[x] > ref_stats[pop]['lcl']]
                vertical= vertical[:1]

                pop_thresh_dict[pop_name][stat_request]= vertical[0]

                ver_bar= [.99,1.002]
                if len(vertical):
                    plt.plot(vertical * 2,ver_bar,color='red',linestyle=style)

                if plot_var:
                    plt.errorbar(surface,y,yerr=error,label= 'var',color= pop_col)   
                else:
                    plt.plot(surface,y,color= pop_col)

                line_surf= [min(counts_dict[pop]['sizes']),max(counts_dict[pop]['sizes'])]

                for z in ['mean','ucl','lcl']:
                    
                    y= [ref_stats[pop][z]] * 2
                    if z == 'mean': 
                        plt.plot(line_surf,y,linestyle=style, 
                            label= '{}'.format(pop_name),color= pop_col)
                    else:
                        plt.plot(line_surf,y,linestyle=style,
                            color= pop_col)

        plt.xlabel(xlab)
        if set_ylim:
            plt.ylim(*ver_bar)
        plt.ylabel(ylab)
        plt.title('grid SSD / sample proportion - control - {}'.format(pop_name), fontsize= 20)

        plt.legend()
        plt.savefig(fig_dir + sims_target + '_{}_gridSSD_ALL_{}_sbsD{}.pdf'.format(stat_request,["",'wSD'][plot_var],['','ZOOM'][int(set_ylim)]),
            format= 'pdf',bbox_inches='tight')
        plt.close()

    ###########################################
    ###########################################
    #### ii. for pop distribution - histograms
    nbins= 50
    opacity= .7

    rates_show= 3

    for pop in pop_diff_dict.keys():
        pop_name= pop_names_dict[pop]
        surface= counts_dict[pop]['sizes']
        steps= np.linspace(0,len(surface)-1,rates_show,dtype= int)
        #
        dtup= [ref_dict[pop]] + [pop_diff_stats[pop][x] for x in steps]
        dtup_len= [len(x) for x in dtup]

        min_len= min(dtup_len)
        #print(dtup)
        sample_dtup= [np.random.choice(x, min_len,replace= False) for x in dtup]
        #

        plt.figure(figsize=(10,10))

        plt.hist(sample_dtup[0],
             color = 'blue', label= 'full population',bins= nbins, alpha= .8)

        for idx in range(rates_show):
            hvec= sample_dtup[idx+1]
            hvec= np.array(hvec)
            hvec= hvec[~np.isnan(hvec)]
            if len(list(set(hvec))) < 3:
                continue
            plt.hist(hvec, alpha= opacity, bins= nbins,
                label= str(surface[steps[idx]]))

        plt.xlabel("MS sd - nbins: {}".format(nbins),fontsize= 20)

        plt.ylabel("frequency",fontsize= 20)
        plt.title(str(pop),fontsize= 20)

        plt.legend()
        plt.savefig(fig_dir + sims_target + '_{}_MSsd_hist_{}_sbsD.png'.format(stat_request,pop_name),bbox_inches='tight')
        plt.close()


    #############################################
    ############### iii. convergence Across populations. 
    ###############
    ylab= 'mean p-val'
    plt.figure(figsize=(10,10))

    for pop in pop_diff_dict.keys():
        pop_name= pop_names_dict[pop]
        pop_col= pop_colors[pop_name]
        x= counts_dict[pop]['sizes']
        y= [np.log(x) for x in average_cdf[pop]]
        
        plt.plot(x,y,label= pop_name,color= pop_col)

    plt.plot([0,sampling[0]],[np.log(pval_threshold/2)]*2,label= str(pval_threshold))

    plt.xlabel(xlab,fontsize=20)

    plt.ylabel(ylab,fontsize= 20)
    plt.ylim(np.log(1e-4),0.1)
    plt.title('grid SSD / sample proportion - control')

    plt.legend()
    plt.savefig(fig_dir + sims_target + '_{}_Variance_Convergence_sbsD.pdf'.format(stat_request),format= 'pdf',bbox_inches='tight')
    plt.close()


####
#### output threshold pvalue
table_dir= 'metadata'
os.makedirs(table_dir, exist_ok=True)
table_dir= table_dir+"/"

table_file= table_dir + sims_target + '_pval{}_table.txt'.format(pval_threshold)

trail= []
for pop in pop_thresh_dict.keys():
    for st,obs in pop_thresh_dict[pop].items():
        trail.append([st,pop,obs])

trail= np.array(trail,dtype= str)
trail= ['\t'.join(x) for x in trail]
trail= '\n '.join(trail)

with open(table_file,'w') as fp:
    fp.write(trail)
