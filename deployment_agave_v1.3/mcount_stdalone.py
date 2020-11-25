



#from tools.SLiM_pipe_tools import mutation_counter_launch
import os
import numpy as np
import itertools as it
import collections

def recursively_default_dict():
    return collections.defaultdict(recursively_default_dict)


import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--sims', type=str,
                    default='sim1')

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

parser.add_argument('--reps', type=int, # replicates per sample size.
                    default=1)

parser.add_argument('--stepup', type=str, # type of analysis - increment / other.
                    default='increment')

parser.add_argument('--simsN', type=int, # Numberr of simulations to read from available. random. 
                    default=0)           # if 0 uses all sims in sims_dir. 

parser.add_argument('--ploidy', type=int, # 
                    default=2)

parser.add_argument('--cwd', type=str,
                    default='./')

parser.add_argument('--db', type=str,
                    default='out_db.txt')

parser.add_argument('--collapsed', action= 'store_true', # collapse mutation counts.
                    default=False)

parser.add_argument('--freqs', action='store_true', # return allele frequencies.
                    default=False)

parser.add_argument('--haps', action='store_true', # analyse by haplotype.
                    default=False)

parser.add_argument('--diffs', action='store_true', # read file of ancestral SNPs and polarise.
                    default=False)


args = parser.parse_args()


#####
#os.chdir(args.cwd)
#####

from tools.input_stdalone import (
	MC_sample_matrix_stdlone
	)

from tools.mcounter_tools import (
    mcounter_deploy_v2
    )

from tools.fasta_utilities import (
    get_mutations, kmer_comp_index, kmer_mut_index
    )

from tools.compare_utilities import (
    gzip_request
)

from INFO_db import INFO_dict

########
########

## directories
main_dir= os.getcwd() + '/'
count_dir= main_dir + 'mutation_counter/count/'
dir_launch= main_dir + 'mutation_counter'
muted_dir= main_dir + 'mutation_counter/data/mutation_count/'
sims_dir= INFO_dict[args.species]['dirs'][args.data]
sims_target= sims_dir.split('/')[-2]
mutlog= 'toMut.log'

indfile= INFO_dict[args.species]['w_indfile']['sims']

pop_names_dict= {g:z for z,g in INFO_dict[args.species]['pop_dict'].items()}

######################################
######################################
##### READ DATA

sampling= [args.samp,args.steps,args.reps]
sampling_str= str(sampling[0])

row= [64,32][int(args.collapsed)]
col= 3
ksize= 3
bases= 'ACGT'

#####
#####
#####

sims_list= args.sims.split(',')

for sim in sims_list:
	print(args.diffs)
	print(args.freqs)

	report= MC_sample_matrix_stdlone(sim,out_db= args.db,min_size= args.minS, samp= sampling, stepup= args.stepup, sim_dir= sims_dir,
                                indfile= indfile, ploidy= args.ploidy,haps_extract=args.haps, ksize= ksize, bases= bases,
	                          diffs= args.diffs,collapsed= args.collapsed,row= row, col= col,
	                       sample_sim= args.simsN,freq_extract= args.freqs)

