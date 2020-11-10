

import os
import numpy as np
import itertools as it

import collections
def recursively_default_dict():
    return collections.defaultdict(recursively_default_dict)

from tools.power_cofactors import (
	sdict_get
	)


import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--species', type=str,
                    default='chimp')

parser.add_argument('--data', type=str,
                    default='sims')

parser.add_argument('--collapsed', action= 'store_true',
                    default= False)

parser.add_argument('--one_tail', action= 'store_false',
                    default= True, help= '')

parser.add_argument('--orelse', action= 'store_true',
                    default= False, help= 'if true, have rate folds in either pop')

parser.add_argument('--norm', action= 'store_true',
                    default= False, help= 'return Z stat cdf')

parser.add_argument('--maxS', type=int,
                    default=50, help= 'max samp range')

parser.add_argument('--Nsteps', type=int,
                    default=0, help= 'if given, subset size steps by Nsteps.')

parser.add_argument('--Rsteps', type=int,
                    default=5, help= '# of synthetic rates in fold rates if no reference in fold_dict')

parser.add_argument('--rangeF', type=str,
                    default='.9,2', help= 'comma delimited range for synthetic rates.')


parser.add_argument('--levels', type=int,
                    default= 6, help= "contour plot levems.")

parser.add_argument('--clim', type=str,
                    default= "", help= "contour range, comma separated float")

parser.add_argument('--width', type=int,
                    default= 12, help= "savgol windows size. must be uneven.")

parser.add_argument('--height', type=int,
                    default= 10, help= "savgol windows size. must be uneven.")

###
### OPTIONAL

parser.add_argument('--savgol', action= 'store_true',
                    default= False, help= 'if given performs savgol filter')

parser.add_argument('--svg_w', type=int,
                    default= 55, help= "savgol windows size. must be uneven.")

parser.add_argument('--meanS', action= 'store_false',
                    default= False, help= "if true, proxy normal nPA sampling")

parser.add_argument('--Nreps', type=int,
                    default= 50, help= "nsamp in case of meanS")

args = parser.parse_args()


from INFO_db import INFO_dict

species= args.species

pop_names= INFO_dict[species]['pop_dict']
pop_names= {g:z for z,g in pop_names.items()}

figdir= 'Figures/{}/'.format(species)


###########################################
###########################################
## parameters lazy
ref_pop= INFO_dict[species]['ref']

# synthetic fold changes argument
rate_bornes= args.rangeF
rate_bornes= rate_bornes.split(",")
rate_bornes= [float(x) for x in rate_bornes]

# sample size arguments for contour plots.
Nsteps= args.Nsteps # subset pop sizes for contour plots (linrange)
sample_range= [2,args.maxS] # min-max
samp_max= sample_range[1]

# analysis arguments.
orelse= args.orelse # if True, fold shifts are assumed to have happened in either pop. if False, only in non-reference pop.
pop_ref_use=  [ref_pop,""][int(orelse)] 

norm_z= args.norm # return z statistic or not.
min_t= [0,.5][int(norm_z)]

two_tail= args.one_tail # one or two tail test.


# fold rate arguments for line plots.
zoom_lib= {
    'extended': {
        'rate_min': 1e-4,
        'rate_max': 8,
        'rate_steps': 1000
    },
    'zoom': {
        'rate_min': .8,
        'rate_max': 3,
        'rate_steps': 1000
    }
}

savgol= args.savgol
svg_w= args.svg_w


# plotting arguments.
height= args.height
width= args.width


#############################
############################# read db_counts
from tools.dbcount import (
    db_collect
)

db_dir= INFO_dict[species]['dirs'][args.data]
db_neigh= os.listdir(db_dir)

pops_dict, count_props_array, muts= db_collect(db_neigh, pop_names,db_dir= db_dir,
												mean_samp= args.meanS, Nreps= args.Nreps)
print('count_props_array shape:')
print(count_props_array.shape)

pop_size_dict= {z: sorted(list(g.keys())) for z,g in pops_dict.items()}

######################
###################### MUTATION PROPORTIONS
### collapse counts in count_props_array (FP counts by sim array).
### only used here to get the distribution of count proportions across simulations / pops. 
### could be used for more. 

if args.collapsed:
	from tools.dbcount import collapse_mat

	count_props_array, kmers= collapse_mat(count_props_array,muts)

###
### median count proportions, unfolded.
count_props_array= np.mean(count_props_array,axis= 0)
prob_median= np.median(count_props_array)


## reference proportion
p= prob_median # reference probability.
print('reference probability:')
print(p)


########################
########################
### fold change dictionary
###
###

from RATE_dict import (
	fold_dict, synth_fold
	)

###
### Power analyses
species_fold= fold_dict[species]

### make a synthetic rate range if none are available.
fold_dict= synth_fold(fold_dict, species, ref_pop, pop_names,
	rate_bornes= rate_bornes, rate_steps= args.Rsteps)



#############################################
############################################# 
####
####
#### CONTOUR
####
from factor_tools.factor_tools import (
    get_sizes
)

from tools.power_cofactors import (
	rate_by_size, comb_power_traverse
	)

from factor_tools.factor_plots import plot_contourPLT

clim= []

if args.clim:
	clim= [float(x) for x in args.clim.split(",")]

for pop_skew,avail_folds in species_fold.items():
	comb_select= tuple([ref_pop,pop_skew])
	if len(set(comb_select)) == 1:
		continue

	avail_folds= {z:g for g,z in avail_folds.items()}
	rate_range= sorted(avail_folds.keys())

	bi_select, samp_order= get_sizes(comb_select, pop_size_dict,Nsteps= Nsteps, samp_max= samp_max)

	rate_stats= rate_by_size(comb_select, pops_dict, rate_range, bi_select, p= p, 
				two_tail= two_tail, orelse= orelse, norm_z= norm_z)

	for rate in rate_range:

		figname= 'Contour_{}_r{}_{}-{}.pdf'.format(species,rate,comb_select[0],comb_select[1])
		figname= figdir + figname

		title_add= ' pops: {} ; ref p= {}'.format(comb_select,p)
		if rate in avail_folds.keys():
			title_add= title_add + '; r= {}'.format(avail_folds[rate])

		xtitle= 'N ' + comb_select[0]
		ytitle= 'N ' + comb_select[1]

		plot_contourPLT(rate_stats,rate,func_select= np.nanmean,
		                 height= height,width= width+2,title_add= title_add, xtitle= xtitle, ytitle= ytitle,
		                 savefig=figname,clim= clim, levels= args.levels)





##################################################################
##################################################################
#####
###
### Power folds.
#####
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline
from scipy.signal import savgol_filter


####
#### Full
zoom_type= 'extended'

for size_sel in ["Biblio", "Fig1"]:

	sdict= sdict_get(INFO_dict, pops_dict, size_sel= size_sel, species= species)
	
	stats_comp, comb_sdict= comb_power_traverse(pops_dict,sdict,p= p, 
	                                two_tail= two_tail, **zoom_lib[zoom_type], 
	                                orelse= orelse, norm_z= norm_z, pop_ref= pop_ref_use)

	figname= 'power_1M_{}_{}_ref.{}.pdf'.format(species,size_sel,zoom_type,ref_pop)
	figname= figdir + figname

	plt.figure(figsize=(width,height))

	for comb in stats_comp.keys():
	    X= comb_sdict[comb]
	    Y= stats_comp[comb]

	    if savgol:
	    	Y = savgol_filter(Y,svg_w,3)

	    label= '-'.join(comb)
	    plt.plot(X,Y,label= label)
	    
	plt.xlabel('fold rate shift')
	plt.ylabel('abs Z')

	plt.title('size_ref: {}; species: {}'.format(size_sel,species))

	plt.xticks(fontsize=10)
	plt.yticks(fontsize=10)

	plt.legend()
	plt.savefig(figname)
	plt.close()




###################
###################
### Focus & rates
###

cols_list= ['g','r','c','m','y']

zoom_type= "zoom"
size_sel= "Biblio"

sdict= sdict_get(INFO_dict, pops_dict,size_sel= size_sel, species= species)

stats_comp, comb_sdict= comb_power_traverse(pops_dict,sdict,p= p, 
                                two_tail= two_tail, **zoom_lib[zoom_type], 
                                orelse= orelse, norm_z= norm_z, pop_ref= pop_ref_use)

skew_power_table= []

margin= 0.01

min_Y= min([min(x) for x in stats_comp.values()]) - margin
max_Y= max([max(x) for x in stats_comp.values()]) + margin

for pop_skew,avail_folds in species_fold.items():

	comb_select= tuple([ref_pop,pop_skew])

	avail_folds= {z:g for g,z in avail_folds.items()}

	comb= [x for x in stats_comp.keys() if len(list(set(x).intersection(comb_select))) == len(comb_select)][0]

	X= comb_sdict[comb]
	Y= stats_comp[comb]

	figname= 'power_1M_{}_{}_{}_ref_{}_skew_{}.pdf'.format(species,size_sel,zoom_type,ref_pop,pop_skew)
	figname= figdir + figname

	plt.figure(figsize=(10,10))

	if savgol:
		Y = savgol_filter(Y,svg_w,3)

	label= '-'.join(comb)
	plt.plot(X,Y,label= label)
	plt.ylim(min_Y,max_Y)

	d= 0

	for rate in avail_folds.keys():
		xrate= np.array(X) - rate
		xrate= np.argmin(abs(xrate))
		yrate= Y[xrate]
		xrate= X[xrate]

		skew_power_table.append(
			[label,rate,yrate]
			)

		label_rate= '{} - {}'.format(avail_folds[rate],round(rate,3))

		plt.plot([xrate,xrate],[min(Y),yrate],linestyle= 'dashed', label= label_rate, color= cols_list[d])
		plt.plot([zoom_lib[zoom_type]["rate_min"],xrate],[yrate,yrate],
			linestyle= 'dashed',color= cols_list[d])
		d += 1

	    
	plt.xlabel('fold rate shift')
	plt.ylabel('abs Z')

	plt.title('size_ref: {}; species: {}'.format(size_sel,species))

	plt.xticks(fontsize=10)
	plt.yticks(fontsize=10)

	plt.legend()
	plt.savefig(figname)
	plt.close()


skew_power_table= np.array(skew_power_table,dtype= str)
filename= figdir + 'power_table_ref_{}.txt'.format(ref_pop)

with open(filename,'w') as fp:
	fp.write('\t'.join(['pop_comp','rate','power']) + '\n')

	fp.write('\n'.join(['\t'.join(x) for x in skew_power_table]))