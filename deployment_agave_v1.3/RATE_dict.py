
import numpy as np

fold_dict= {
    "human": {
        'EUR': {
            'TCC>T': 1.56,
            'ACC>T': 1.2,
            'TCT>T': 1.17,
            'CCC>T': 1.06
        },
        'ASN': {
            'TCC>T': 1.00,
            'ACC>T': 0.93,
            'TCT>T': 0.98,
            'CCC>T': 0.95
        }
    },
    'chimp': {}
}



def synth_fold(fold_dict, species, ref_pop, pop_names,
    rate_bornes= [.8,2], rate_steps= 5):
    '''
    synthesise rates between against ref species if empty species dict in fold_dict.
    '''

    if not fold_dict[species]:
        for pop,pop_c in pop_names.items():
            
            if pop_c != ref_pop:
                rate_range= np.linspace(rate_bornes[0],rate_bornes[1],rate_steps,dtype= float)

                fold_dict[species][pop_c]= {
                    "synth"+str(idx): rate_range[idx] for idx in range(len(rate_range))
                }

    return fold_dict
