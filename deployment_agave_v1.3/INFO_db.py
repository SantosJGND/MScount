
INFO_dict= {
    'test': {
        'dirs': {
            'sims': 'data/'
        },
        'w_indfile': {
            'rd': 'chimp_ind_assignments.txt',
            'sims': 'ind_assignments.txt'
        },
        'g_indfile': '/home/jgarc235/Chimp/chimp_ind_assignments.txt',
        'pop_dict': {
            'schweinfurthii': 'pop0',
            'troglodytes': 'pop1',
            'ellioti': 'pop2',
            'verus': 'pop3'
        },
        'pop_colors': {
            'schweinfurthii': 'red',
            'troglodytes': 'blue',
            'ellioti': 'green',
            'verus': 'orange'
        },
        'assembly': 'panTro5'
    }, 
    'chimp':{
        'ref': 'verus',
        'dirs':{
            'rd': '/home/jgarc235/Chimp/vcf_data/chimp_1MB/',
            'FP': '/scratch/jgarc235/SLiM/chimp/FP/', 
            '1K': '/home/jgarc235/SLiM/PM13_chimps/mutation_counter/data/pm13_4A_1m/',
            '50M': '/scratch/jgarc235/SLiM/chimp/50M/',
            'local_1K': 'D:/GitHub/fine-scale-mutation-spectrum-master/slim_pipe/mutation_counter/data/pm13_4A_1m/',
            'local_counts': 'D:/GitHub/SLiM/Analyses/Local_joao/db_dir/pm13_4A_1m__population/'
        },
        'w_indfile': {
            'rd': 'chimp_ind_assignments.txt',
            'sims': 'ind_assignments.txt'
        },
        #'g_indfile': '/home/jgarc235/Chimp/chimp_ind_assignments.txt',
        'g_indfile': 'D:/GitHub/SLiM/Analyses/Local_joao/metadata/chimp_ind_assignments.txt',
        'pop_dict': {
            'schweinfurthii': 'pop0',
            'troglodytes': 'pop1',
            'ellioti': 'pop2',
            'verus': 'pop3'
        },
        'pop_colors': {
            'schweinfurthii': 'red',
            'troglodytes': 'blue',
            'ellioti': 'green',
            'verus': 'orange'
        },
        'pop_style': {
            'schweinfurthii': ':',
            'troglodytes': '-.',
            'ellioti': '--',
            'verus': '-'
        },
        'assembly': 'panTro5'
    },
    'human':{
        'ref': 'AFR',
        'dirs':{
            'rd': '/home/jgarc235/Human/sim_compare/data/phase1_1MB/',
            'FP': '/scratch/jgarc235/SLiM/demos/data/gravel_FP/',
            '1K': '/home/jgarc235/SLiM/demos/mutation_counter/data/gravel1_1MB/',
            '50M': '/scratch/jgarc235/SLiM/demos/data/50M/',
            'local_1K': 'D:/GitHub/fine-scale-mutation-spectrum-master/slim_pipe/mutation_counter/data/gravel_1m_1Ksamp/',
            'local_counts': 'D:/GitHub/SLiM/Analyses/Local_joao/db_dir/gravel_1m_1Ksamp__population/'
            },
        'w_indfile': {
            'rd': 'integrated_call_samples.20101123.ALL.panel_regions.txt',
            'sims': 'ind_assignments.txt'
            },
        #'g_indfile': '/home/jgarc235/Human/ind_assignment/integrated_call_samples.20101123.ALL.panel_regions.txt',
        'g_indfile': 'D:/GitHub/SLiM/Analyses/Local_joao/metadata/integrated_call_samples.20101123.ALL.panel_regions.txt',
        'pop_dict': {
            'AFR': 'pop0',
            'EUR': 'pop1',
            'ASN': 'pop2'
        },
        'pop_colors': {
            'AFR': 'blue',
            'EUR': 'red',
            'ASN': 'green'
        },
        'pop_style': {
            'AFR': ':',
            'EUR': '-.',
            'ASN': '--'
        },
        'assembly': 'hg38'
    },
    'rhesus':{
        'dirs':{
            'rd': '/home/jgarc235/Rhesus/Windows/rhe18_1M/',
            'sims': '/home/jgarc235/SLiM/Rhesus_Liu18/mutation_counter/data/rheLiu18_1M/'
            },
        'w_indfile': {
            'rd': 'ind_assignments.txt',
            'sims': 'ind_assignments.txt'
            },
        'g_indfile': '/home/jgarc235/Rhesus/Mutation_study/ind_assignments.txt',
        'pop_dict':{
            'tcheliensis': 'pop0',
            'littoralis': 'pop1',
            'brevicaudus': 'pop2',
            'lasiotis': 'pop3',
            'mulatta': 'pop4'
        },
        'pop_colors':{
            'tcheliensis': 'blue',
            'littoralis': 'red',
            'brevicaudus': 'green',
            'lasiotis': 'orange',
            'mulatta': 'blueviolet'
        },
        'pop_colors':{
            'tcheliensis': ':',
            'littoralis': '-.',
            'brevicaudus': '--',
            'lasiotis': '-',
            'mulatta': 'densely dotted'
        },
        'assembly': 'rheMac10'
    }
}

