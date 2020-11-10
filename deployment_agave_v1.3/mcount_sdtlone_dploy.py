#####
from tools.input_stdalone import (
	deploy_countDB
	)

from tools.fasta_utilities import (
    get_mutations, kmer_comp_index, kmer_mut_index
    )

from tools.compare_utilities import (
    gzip_request
)


#from tools.SLiM_pipe_tools import mutation_counter_launch
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

parser.add_argument('--reps', type=int, # replicates per sample size.
                    default=1)

parser.add_argument('--stepup', type=str, # type of analysis - increment / other.
                    default='increment')

parser.add_argument('--ploidy', type=int, # 
                    default=1)

parser.add_argument('--mem', type=str,
                    default='15GB')

parser.add_argument('-t', type=str,
                    default='30:00:00')

parser.add_argument('--nodes', type=int,
                    default= 4)

parser.add_argument('--npacks', type=int,
                    default= 1)

parser.add_argument('--simsN', type=int, # Numberr of simulations to read from available. random. 
                    default=0)           # if 0 uses all sims in sims_dir. 

parser.add_argument('--collapsed', action='store_true', # collapse mutation counts.
                    default=False)

parser.add_argument('--haps', action='store_true', # analyse by haplotype.
                    default=False)

parser.add_argument('--debug', action= 'store_true',
                    default= False)

parser.add_argument('--logdir', type=str,
                    default='./log')

parser.add_argument('--deployment', type=str,
                    default='local')


args = parser.parse_args()


#####
if args.simsN > 0:
    args.simsN += 1
#####
from INFO_db import INFO_dict

### should become arguments
ksize= 3
bases= 'ACGT'

######## ###### could be imported from library. 
######## ###### if mcount_stdalone.py diverges.
daughter_script= os.getcwd() + '/mcount_stdalone.py'


sims_dir= INFO_dict[args.species]['dirs'][args.data]
sims_target= sims_dir.split('/')[-2]

request_processed= gzip_request(dir_check= sims_dir,str_format= '',requested= ['.vcf','fa'], func= 'gzip')

#####
##### Prepare database:
db_dir= args.db_dir
sub_dir= '/{}_{}_{}'.format(sims_target,['','coll'][int(args.collapsed)],args.stepup)
db_dir= db_dir + sub_dir
os.makedirs(db_dir, exist_ok=True)
db_dir= db_dir + '/'

db_name= sims_target + '_{}' + '_{}_{}_countDB.txt'.format(['','coll'][int(args.collapsed)],args.stepup)
db_file= db_dir + db_name

####
#### prepare deployment by simulation batch. 
### remove unnecessary arguments. prepare to pass to mcount_stdalone.
stdalone_absent= ['db_dir','npacks','t','mem','nodes',"debug",
                    "deployment", "logdir"]

elements= vars(args)
elements= {z:g for z,g in elements.items() if z not in stdalone_absent}

print(elements)

### command to be passed to batch. 
command_base= 'python -u {} '.format(daughter_script)
command_dir= daughter_script.split('/')[:-1]
command_dir= '/'.join(command_dir) + '/'

###
###
deploy_countDB(elements, command_base, npacks= args.npacks,sim_dir= sims_dir,
                mem= args.mem, t= args.t, nodes= args.nodes,sample_sim= args.simsN,
                out_db= db_file,command_dir= command_dir,debug= args.debug, 
                deployment= args.deployment, log_dir= args.logdir)


