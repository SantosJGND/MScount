import numpy as np
import os
#############################
#############################
from factor_tools.factor_tools import (
    factor_process
)



##################

def read_simCounts(simdb,tag_split= 'C',pop_tag= '_ss'):
	'''
	read file of individual mutation type counts. 
	first 3 columns= simID, pop, ind. 
	header= True.
	'''
	print(simdb)
	with open(simdb,'r') as fp:
		counts= fp.readlines()

	header= counts[0]
	muts= header.strip().split('\t')[3:]
	counts= [x.strip().split('\t') for x in counts[1:]]
	counts= np.array(counts)
	pop_names= counts[:,1]
	for idx in range(len(pop_names)):
		pop= pop_names[idx]
		if pop_tag in pop:
			pop= pop[len(pop_tag):].split('.')[0]
			pop_names[idx]= pop 

	counts[:,1]= pop_names
	
	return counts, muts, header




def info_array_collect(db_dir, 
					row= 24,col= 4, tag_split= 'C', tag_sim= '_ss'):
	'''
	extract mutation count sub-samples from across simulation count dbs. 
	'''

	list_neigh= os.listdir(db_dir)
	lines= []

	for dbf in list_neigh[:1]:
		simdb= db_dir + dbf
		print('reading {}'.format(simdb))
		counts, muts, header= read_simCounts(simdb)

		batch= counts[0][0]
		print('base: {}'.format(batch))
		print(counts.shape)

		lines.append(counts)

	lines= np.concatenate(tuple(lines),axis= 0)

	info_array= lines[:,:3]

	counts= lines[:,3:]
	counts= np.array(counts,dtype= int)

	return info_array, muts, counts




def sim_countsSize(simdb,fact_funcs= [], fact_args= [], fact_names= ['pop','size']):
	'''
	read simulation specific count array, return dictionary of counts per samp size per population. 
	'''

	counts, muts, header= read_simCounts(simdb,tag_split= 'C',pop_tag= '_ss')

	info_array= counts[:,:3]
	counts= counts[:,3:]
	counts= np.array(counts,dtype= int)

	###
	counts_full= counts[:4]
	count_prox= np.array(counts_full)
	count_prox[count_prox == 0]= 1
	
	count_props= counts_full.T / np.sum(count_prox, axis= 1)
	##
	## factor plot 
	# functions
	
	pop_array= info_array[:,1:]

	pop_array, diffs, pops_dict= factor_process(pop_array, np.array(counts), fact_names= fact_names, 
	                                        fact_functs= fact_funcs, fact_args= fact_args)


	return pops_dict, count_props, muts

###############################
###############################

from factor_tools.factor_tools import (
    empty_func, varfac_collect, Nseg_counts
)


def db_collect(db_neigh, pop_names, db_dir= './',
				mean_samp= False, Nreps= 50):
	'''
	Getting segregating PA counts by sample size and population. 
	pops_dict= list of NPA by pop:{samp_size:[]};
	count_props= np.array of FPs across simulations;
	muts= header of the last sim surveyed.

	pops_dict= sizes by samp and pop dictionary.
	If samp_mean: samples from mean and var of nPA list by pop/samp_size assuming normality.
	i.e. if samp counts change per pop/samp size in database. 
	'''

	fact_funcs= [empty_func,varfac_collect]

	fact_names= ['pop','size']

	fact_args= [{},{'func_plot': Nseg_counts}]

	###
	###

	pops_dict= {z: {} for z in pop_names.values()}
	count_props_array= []

	for idx in range(len(db_neigh)):
		simdb= db_dir + db_neigh[idx]

		pops_counts, count_props, muts= sim_countsSize(simdb,fact_funcs= fact_funcs, fact_args= fact_args, fact_names= fact_names)
		pops_counts= {pop_names[z]:g for z,g in pops_counts.items()}

		count_props_array.extend(count_props.T)

		for pop,g in pops_counts.items():
			for si,cts in g.items():
				if si not in pops_dict[pop].keys():
					pops_dict[pop][si]= list(cts)
				else:
					pops_dict[pop][si].extend(cts)

	if mean_samp:

		pops_dict= {
			pop: {
				si: {
					'mean': np.median(t),
					'std': np.std(t)
				} for si,t in pops_dict[pop].items()
			} for pop in pops_dict.keys()
		}

		pops_dict= {
			pop: {
				si: norm.rvs(loc= t['mean'],scale= t['std'], size= Nreps) for si,t in pops_dict[pop].items()
			} for pop in pops_dict.keys()
		}

	count_props_array= np.array(count_props_array)

	return pops_dict, count_props_array, muts


#####
##### COLLAPSE DBCOUNT MATRIX

from tools.fasta_utilities import (
	kmer_comp_index
	)


def collapse_mat(count_props_array,muts):
	'''
	collapse counts_props_array.
	'''
	muts_tuple= [x.split('_') for x in muts]
	mutations= [tuple(x) for x in muts_tuple]
	kmers, kmer_idx= kmer_comp_index(mutations)

	muts_idx= [kmers[''.join(x)] for x in mutations]

	ncount_array= np.zeros((count_props_array.shape[0],len(kmer_idx)))

	for idx in range(len(muts_idx)):
		nidx= muts_idx[idx]
		
		ncount_array[:,nidx]+= count_props_array[:,idx]

	count_props_array= ncount_array

	return count_props_array, kmers