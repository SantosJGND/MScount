##
# Deployment functions.
##
import time
import numpy as np
import os

from tools.SLiM_pipe_tools import (
	sbatch_launch, bash_launch
	)

from tools.input_utilities import (
	VCF_read_filter
	)

from tools.input_cofactors import (
	count_popKmers, process_dir, get_pop_dict,
	ind_assignment_scatter_v1, parse_PA
	)

from tools.fasta_utilities import (
    get_mutations, kmer_comp_index, kmer_mut_index
    )



def get_bins(array_len= 10, npacks= 1):
    '''
    get bins from array.
    '''
    if npacks > array_len:
    	npakcs = array_len

    t_range= np.arange(0,array_len,npacks,dtype= int)

    bins= []
    for idx in range(len(t_range)):
        if idx == len(t_range)-1:
            if t_range[idx] == array_len-1:
                bins[-1][1] = array_len-1
            else:
                bins.append([t_range[idx],array_len-1])
        else:
            bins.append([t_range[idx],t_range[idx+1]-1])

    return bins


def deploy_countDB(elements, command_base, npacks= 1, sim_dir= './',
					mem= '3GB', t= '00:30:00', nodes= 2, sample_sim= 0,
					out_db= 'out_db',command_dir= './',debug= False, deployment= 'local',
					log_dir= './log'):
	"""
	deployment of smaller subpops. 
	- elements: dict, script arguments to pass to daugher script.
	- adds arguments sims, db and cwd to arguments.
	- deploy mcount_standalone on a a number npacks of simulations at a time. 
	- use sbatch_launch() to write bash file to launch mcount_stdalone.py.
	- mcount_stdalone and mcount_stdalone_deploy both read INFO_db for the simulation directory 
	using args.species, passed here in elements arg dict.
	"""
	ti= time.time()
	sims= process_dir(sims_dir= sim_dir)
	print('available {}'.format(len(sims)))

	data_kmer= {}
	data_freqs= {}

	if sample_sim == 0:
	    sample_sim= len(sims)

	print('sample {}'.format(sample_sim))
	sim_sub= np.random.choice(sims,sample_sim,replace= False)


	bins= get_bins(array_len= len(sim_sub),npacks=npacks)
	for cycle in bins:
		
		if len(list(set(cycle))) == 1:
			cycle[1] += 1

		aprehend= sim_sub[cycle[0]:cycle[1]]
		bash_name= '_'.join(aprehend)
		#
		aprehend= ','.join(aprehend)
		arg_dict= dict(elements)

		command_here= str(command_base)

		arg_dict['sims']= aprehend
		arg_dict['db']= out_db
		arg_dict['cwd']= command_dir 

		for arg,val in arg_dict.items():
			suff= '--'
			if len(arg) == 1:
				suff = '-'

			if isinstance(val,bool):
				if val:
					command_here += suff + arg + ' '
			else:
				command_here += suff + arg + ' {} '.format(str(val))

		if deployment == 'local':
			bash_launch(command_here,bash_name,log_dir= log_dir)

		else:
			sbatch_file= sbatch_launch(command_here,bash_name,batch_dir= '', 
				mem= mem,t= t,nodes= nodes,debug= debug)

			os.system('sbatch ' + sbatch_file)


##
## STANDALONE FUNCTIONS
##
##

def MC_sample_matrix_stdlone(sim,out_db= 'out_db.txt',min_size= 80, samp= [5,20,10], stepup= "increment", diffs= False, frequency_range= [0,1],indfile= 'ind_assignments.txt', 
					outemp= 'ind_assignments{}.txt',chrom_idx= 0, prop_gen_used= 1, 
	                sim_dir= 'mutation_counter/data/sims/', segregating= False, scale_genSize= False,
	                outlog= 'indy.log', row= 24,col= 4, single= True, exclude= False, print_summ= False, sample_sim= 0,
	                collapsed= True,bases= 'ACGT',ksize= 3,ploidy= 2, freq_extract= False, sim_del= 'C', tag_sim= '_ss',
	                genome_size= 1,haps_extract= False, return_private= True):
	'''

	'''

	### prepare sim_specific db
	###
	### mutation labels
	mutations= get_mutations(bases= bases,ksize= ksize)
	kmers, kmer_idx= kmer_comp_index(mutations)

	mut_lib= kmer_mut_index(mutations)
	if collapsed:
	    labels= [kmer_idx[x][0] for x in sorted(kmer_idx.keys())]

	else:
	    labels= ['_'.join(x) for x in mutations]

	# write:
	out_db= out_db.format(sim)
	db_dir= out_db.split('/')
	db_name= db_dir[-1]
	db_dir= '/'.join(db_dir[:-1])+'/'
	db_neigh= os.listdir(db_dir)

	header= ['SIM','POP','N'] + labels
	header= '\t'.join(header) + '\n'

	if db_name not in db_neigh:
	    with open(out_db,'w') as fp:
	        fp.write(header)
	###
	###
	ti= time.time()

	## chromosome
	chrom= sim.split('.')[chrom_idx].split(sim_del)[-1].strip('chr')
	pop_dict, inds= get_pop_dict(sim,dir_sim= sim_dir,indfile= indfile,haps_extract= haps_extract,
		return_inds= True)
	if haps_extract:
		inds= list(inds) + list(inds)
		inds= np.array(inds)

	print(len(inds))
	print(inds[:10])

	total_inds= sum([len(x) for x in pop_dict.values()])

	if exclude:
	    files= read_exclude()
	else:
	    files= {}
	
	data_kmer= {}
	tags= []
	sim_extend= []
	chroms= []


	### read vcf
	t0= time.time()
	Window, mut_matrix, scale= VCF_read_filter(sim, sim_dir= sim_dir,chrom= chrom,haps_extract= haps_extract, scale_genSize= scale_genSize,
	    collapsed= collapsed,min_size= min_size, samp= samp, stepup= stepup, outemp= outemp,
	    indfile= indfile,diffs= diffs,bases= bases, ksize= ksize, ploidy= ploidy)

	mut_idx= np.argmax(mut_matrix,axis= 0)

	tag_list, tag_dict, pop_dict= ind_assignment_scatter_v1(sim,dir_sim= sim_dir,
	                  min_size= min_size, samp= samp, stepup= stepup, outemp= outemp,indfile= indfile,
	                  inds= inds, pop_dict= pop_dict,pop_sub= True)

	t1= time.time()
	read_time= t1- t0
	if not len(Window) or Window.shape[0] < total_inds:
	    return ''

	## counts for no tag sim:
	s0= time.time()
	PA_dict= parse_PA(Window, pop_dict,frequency_range= frequency_range)

	t2= time.time()
	print('time elapsed ref: {} m'.format((t2 - t1)/60))


	new_wind= {}
	for pop in pop_dict.keys():

		klist= list(pop_dict[pop])
		private= list(PA_dict[pop])
		pop_wind= Window[klist,:]
		pop_wind= pop_wind[:,private]

		##
		new_wind[pop]= {
			'array': pop_wind,
			'mut': [mut_idx[x] for x in private]
			}

		##
		#

	for pop in pop_dict.keys():
		dict_pop= {pop: list(range(len(pop_dict[pop])))}
		pop_summary, dummy= count_popKmers(new_wind[pop]['array'], mut_matrix, new_wind[pop]['mut'], dict_pop, 
												row=row,col=col)

		nline= [sim,pop,len(pop_dict[pop])]
		ncounts= pop_summary['counts'][pop]
		ncounts= ncounts.reshape(1,np.prod(ncounts.shape))[0]
		nline= nline + list(ncounts)
		nline= [str(x) for x in nline]
		with open(out_db,'a') as fp:
			fp.write('\t'.join(nline) + '\n')

	t3= time.time()
	print('time elapsed ref PA: {} m'.format((t3 - t2)/60))

	if len(tag_list):
		###
		print('len tag list: {}'.format(len(tag_list)))
		###

		for idx in range(len(tag_list)):
			
			tag= tag_list[idx]

			pop_dict= tag_dict[tag]
			pop_ori= list(pop_dict.keys())[0]
			if tag_sim in pop_ori:
			    pop_ori= pop_ori[len(tag_sim):].split('.')[0]

			pop_summary, dummy= count_popKmers(new_wind[pop_ori]['array'], mut_matrix, new_wind[pop_ori]['mut'], pop_dict,
						row=row,col=col,segregating= True,single= True)

			for pop in pop_summary['counts'].keys():
				nline= [sim,pop,len(pop_dict[pop])]
				ncounts= pop_summary['counts'][pop]
				ncounts= ncounts.reshape(1,np.prod(ncounts.shape))[0]
				nline= nline + list(ncounts)
				nline= [str(x) for x in nline]
				with open(out_db,'a') as fp:
					fp.write('\t'.join(nline) + '\n')

	return ''



