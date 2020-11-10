
import os
import numpy as np 
from scipy import stats

import itertools as it

import collections
def recursively_default_dict():
    return collections.defaultdict(recursively_default_dict)


def sdict_get(INFO_dict, pops_dict, size_sel= 'Fig1', species= 'human'):
	'''
	get fixed population sizes, from either reference or a ind-pop file.
	'''

	if size_sel== 'Fig1':
	    keyget= 'ref'
	    sdict= sconv_get(species)
	    sdict= sdict[keyget]

	else:
	    ind_file= INFO_dict[species]['g_indfile']
	    sdict= get_pop_dict(indfile= ind_file)
	    sdict= {z:len(g) for z,g in sdict.items() if z in pops_dict.keys()}

	return sdict


def sconv_get(species,metadir= './metadata', key= "FP",pval= 0.01):
    
    dr= recursively_default_dict()
    
    filename= metadir+ '/{}_{}_pval{}_table.txt'.format(species,key,pval)
    
    with open(filename,'r') as fp:
        lines= fp.readlines()
    
    lines= [x.strip().split() for x in lines]
    
    for lin in lines:
        dr[lin[0]][lin[1]]= int(lin[2])
    
    return dr

def get_pop_dict(indfile= 'ind_assignments.txt'):

    #ind_assignments= dir_sim + reference + '/' + indfile
    
    with open(indfile,'r') as f:
        inds= f.readlines()
    
    inds= [x.split() for x in inds]
    inds= [x for x in inds if x]
    pops= np.array(inds)[:,1]
    pop_dict= {
        z: [x for x in range(len(pops)) if pops[x] == z] for z in list(set(pops))
    }

    return pop_dict


###
###

def z_calc(p1, p2, n1, n2):
    '''
    calculate Z score.
    n1 and n2 are numpy vectors.
    sets values to 0 without returning error. 
    '''

    sumN= (n1+n2)

    null_comb= (n1 == 0) & (n2 == 0)
    null1= np.array(n1) 
    null1= null1 == 0
    n1[null1]= 1
    null2= np.array(n2) 
    null2= null2 == 0
    n2[null2]= 1

    p_star = (p1*n1 + p2*n2) / (n1 + n2)
    nl= p2 - p1
    
    tl= p_star*(1 - p_star)*((1.0 / n1) + (1.0 / n2))
    tl[null_comb]= 0

    return nl / np.sqrt(tl)

def power_calc(p1,p2,N1,N2, norm= True,two_tail= True):
    '''
    prepare z_calc for the size vectors N1 and N2 given proportions p1 and p2. vectorize.
    '''

    if len(N1) != len(N2):

        min_len= min([len(y) for y in [N1,N2]])
        N1,N2= [np.random.choice(x,min_len,replace= False) for x in [N1,N2]]

    p_array= np.repeat(p1,len(N1))
    pn_array= np.repeat(p2,len(N1))
    
    z= z_calc(p_array, pn_array, n1=N1, n2=N2)

    if two_tail:
        z= abs(z)

    if norm:
        z = stats.norm.cdf(z)

    return z

###
###

def comb_power_traverse(pops_dict,sdict,p= 1/96, rate_min= 1e-4,rate_max= 15, rate_steps= 1000,
							two_tail= True, orelse= True, norm_z= True, pop_ref= ''):

	rate_range= np.linspace(rate_min,rate_max,rate_steps,dtype= float)

	if pop_ref:
		comb_pops= [(pop_ref,x) for x in pops_dict.keys() if x != pop_ref]
	else:
		comb_pops= it.combinations(pops_dict.keys(),2)
		comb_pops= list(comb_pops)

	comb_stats_dict= {}

	for comb_sel in comb_pops:
	    bi_sel= [sdict[x] for x in comb_sel]
	    bi_sel= [tuple(bi_sel)]

	    comb_stats= rate_by_size(comb_sel, pops_dict, rate_range, bi_sel, p= p, 
	                 orelse= orelse, norm_z= norm_z,two_tail= two_tail)

	    comb_stats_dict[comb_sel]= comb_stats[bi_sel[0]]


	stats_comp= {x: [comb_stats_dict[x][z]['pval'][0] for z in sorted(g.keys())] for x,g in comb_stats_dict.items()}
	stats_comp= {z:np.array(g) for z,g in stats_comp.items()}

	sizes_dict= {z: sorted(comb_stats_dict[z].keys()) for z in comb_stats_dict.keys()}

	return stats_comp, sizes_dict



####
####
####

def comb_power(comb_select, pops_dict, ref_idx= 0, bi= 5, Nreps= 20,p1= 1/192,p2= 1/192,
                norm_z= True, two_tail= True):
    '''
    return power for a range of count vectors for a combination of populations comb_select.
    '''
    pop_idx= {
        comb_select[idx]: idx for idx in range(len(comb_select))
    }
    
    pop_skew= comb_select[ref_idx]
    pop_ref= [x for x in comb_select if x != pop_skew][0]

    Nskew= pops_dict[pop_skew][bi[pop_idx[pop_skew]]]
    Nskew= np.array(Nskew)
    #
    Nref= pops_dict[pop_ref][bi[pop_idx[pop_ref]]]
    Nref= np.array(Nref)
    #

    probH2= power_calc(p1,p2,Nref,Nskew, norm= norm_z, two_tail= two_tail)
    
    return probH2


##


def rate_by_size(comb_select, pops_dict,rate_range, bi_select, p= 1/192, 
                 Nreps= 50, two_tail= False, orelse= True, norm_z= True):
    
    '''
    Run analysis using pops_dict stats combinning two populations (comb_select tuple), stats dict has two levels
    Second level values have 'mean' and 'std' values. 
    '''

    rate_stats= recursively_default_dict()

    for bi in bi_select:
        for rate in rate_range:

            new_rate= p * rate

            rate_stats[bi][rate]= {
                x: [] for x in [*comb_select,'pval']
            }

            pops_probs= []

            sel_idx= list(range(len(comb_select)))

            if not orelse:
                sel_idx= [sel_idx[-1]]

            for pop_sel in sel_idx:
                
                probH2= comb_power(comb_select, pops_dict, p1= p, p2= new_rate, two_tail= two_tail,
                                   ref_idx= pop_sel, bi= bi, Nreps= Nreps, norm_z= norm_z)
                
                pops_probs.append(probH2)

            pops_probs= np.array(pops_probs)
            pops_probs= np.nanmean(pops_probs,axis= 1)

            if len(pops_probs):
                for idx in range(len(comb_select)):
                    rate_stats[bi][rate][comb_select[idx]].append(bi[idx])

                rate_stats[bi][rate]['pval'].append(np.min(pops_probs))
    
    return rate_stats








################################################################
################################################################
########
######## Deprecated

def rate_run(comb_select, pops_dict, rate_range, bi_select, p= 1/192, Nreps= 20, alpha= .2,
                orelse= True, norm_z= True):

    rate_min= []

    for bi in bi_select:

        idx= 0
        min_found= -1
        
        while min_found < 0:
            rate= rate_range[idx]
            new_rate= p * rate

            pops_probs= []

            sel_idx= list(range(len(comb_select)))

            if not orelse:
                sel_idx= [sel_idx[-1]]

            for pop_sel in sel_idx:

                probH2= comb_power(comb_select, pops_dict, p1= p, p2= new_rate,
                                   ref_idx= pop_sel, bi= bi, Nreps= Nreps, norm_z= norm_z)

                pops_probs.append(probH2)

            pops_probs= np.array(pops_probs)
            #
            pops_probs= np.nanmean(pops_probs,axis= 1)

            if 1 - np.min(pops_probs) <= alpha / 2 or np.min(pops_probs) <= alpha / 2:
                min_found= rate

            idx += 1

            if idx == len(rate_range):
                min_found= rate
            else:
                rate= rate_range[idx]

        rate_min.append(min_found)
    
    return rate_min


##
##

def comb_FixedPower(comb_pops, pops_dict, p= 1/192, rate_range= [], Nsteps= 20, samp_max= 50,
                 set_unique= True, Nreps= 20, alpha= .05, norm_z= True):
    '''
    run rate_run on combinations of pops given pops_sizes dict.
    '''

    pop_size_dict= {z: sorted(list(g.keys())) for z,g in pops_dict.items()}
    
    if not len(rate_range):
        rate_range= np.linspace(1,20,20)
        
    comb_dict={}

    for comb_select in comb_pops:

        pop_idx= {
            comb_select[idx]: idx for idx in range(len(comb_select))
        }

        bi_select, samp_order= get_sizes(comb_select, pop_size_dict,Nsteps= Nsteps, samp_max= samp_max,
                 set_unique= set_unique)

        rate_min= rate_run(comb_select, pops_dict, rate_range, bi_select,  p= p, Nreps= Nreps, 
        alpha= alpha, norm_z= norm_z)

        comb_dict[comb_select]= {
            "sizes": [x[0] for x in bi_select],
            "rates": rate_min

        }
    
    return comb_dict
