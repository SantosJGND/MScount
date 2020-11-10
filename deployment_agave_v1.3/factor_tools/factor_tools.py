
import numpy as np
import itertools as it
import collections
def recursively_default_dict():
    return collections.defaultdict(recursively_default_dict)


import plotly
#import chart_studio.plotly as py
import plotly.graph_objs as go


##
from scipy.stats import norm
from scipy.stats import binom
import math
import pandas as pd
from scipy import stats


###
###
## math functions
from scipy import stats


###
###


def get_sizes(comb_id, pop_sizes,Nsteps= 0, samp_max= 0,
             set_unique= False):
    '''
    given tuple IDs, {pop: [sizes]} dictionary, create list of size combinations where sizes used are subset:
    - linearly, Nsteps;
    - curtailed below samp_max. <=.
    '''
    sizes= [list(pop_sizes[pop]) for pop in comb_id]
    if samp_max:
        sizes= [[x for x in y if x <= samp_max] for y in sizes]
    
    if Nsteps:
        sizes_subsel= [np.linspace(0,len(x)-1,Nsteps, dtype= int) for x in sizes]
        sizes= [[sizes[y][x] for x in sizes_subsel[y]] for y in range(len(comb_id))]
        sizes= [list(set(x)) for x in sizes]
    
    if set_unique:
        samp_order= set(sizes[0]) & set(sizes[1])
        samp_order= sorted(list(samp_order))
        bi_select= [(x,x) for x in samp_order]
    
    else:
        samp_order= [list(sorted(x)) for x in sizes]
        bi_select= it.product(*samp_order)
        bi_select= list(bi_select)
    
    return bi_select, samp_order

##

def dict_from_array(arrayc):
    '''
    create index dictionary fromm array w/ 1 loop.
    '''
    array_dict= recursively_default_dict()
    for idx in range(len(arrayc)):
        array_dict[arrayc[idx]][idx]= 0
    
    array_dict= {z: sorted(list(g.keys())) for z,g in array_dict.items()}
    
    return array_dict
    

### layer info

# layer functions.
def empty_func(info_array,idx, counts, goods= {}):
    '''
    '''
    return info_array, counts, goods



###
###

def factor_process(info_array, counts, fact_names= [], fact_functs= [], fact_args= [],goods= {}):
    '''
    '''
    for idx in range(info_array.shape[1]):
        
        print('layer: {}'.format(fact_names[idx]))
        func_here= fact_functs[idx]
        info_array, counts, goods= func_here(info_array, idx, counts, goods= goods,**fact_args[idx])
    
    return info_array, counts, goods



####
####
#### INNER FUNCTIONS



def Nseg_counts(array_count, ave= True):
    
    t= np.sum(array_count,axis= 1)
    
    return t
 

def Nseg_dict(array_count):
    
    t= np.sum(array_count,axis= 1)
    
    return {
        'mean': np.mean(t),
        'std': np.std(t)
    }


###
###
### FACTOR FUNCTIONS


def size_la_diffs(info_array,idx,counts,goods= {}):
    '''
    create array of differences to full population proportions corrected for simulation. 
    '''
    sim_dict= info_array[:,(idx-2)]
    sim_dict= dict_from_array(sim_dict)
    
    for sim,sidx in sim_dict.items():
        #
        pop_dict= info_array[sidx,(idx-1)]
        pop_dict= dict_from_array(pop_dict)
        pop_dict= {z: [sidx[x] for x in g] for z,g in pop_dict.items()}
        
        array_dict= {z: counts[g,:] for z,g in pop_dict.items()}

        for pop,pidx in pop_dict.items():
            #
            size_dict= info_array[pidx,idx]
            size_dict= np.array(size_dict,dtype= int)
            size_dict= dict_from_array(size_dict)
            ref= max(size_dict.keys())
            
            ref_diffs= np.mean(array_dict[pop][size_dict[ref],:],axis= 0)
            array_dict[pop]= array_dict[pop] - ref_diffs
        
            counts[pidx]= array_dict[pop]    
    
    
    return info_array,counts, goods






def varfac_plotly(info_array,idx,counts,goods= {},exc_max= True,func_plot= np.nanmean,plot_fig= False):
    '''
    factor function with optional plotly output and dictionary output of INNER FUNCTION, see Nseg_dict.
    '''
    #
    pop_dict= info_array[:,(idx-1)]
    pop_dict= dict_from_array(pop_dict)

    array_dict= {z: counts[g,:] for z,g in pop_dict.items()}
    
    fig= []
    
    
    for pop,pidx in pop_dict.items():
        #
        size_dict= info_array[pidx,idx]
        size_dict= np.array(size_dict,dtype= int)
        size_dict= dict_from_array(size_dict)
        
        sizes= sorted(size_dict.keys())[:-1]
        
        props= [func_plot(array_dict[pop][size_dict[z],:]) for z in sizes]
        prop_means= [x['mean'] for x in props]
        prop_std= [x['std'] for x in props]
        fig.append(go.Scatter(
            x= sizes,
            y= prop_means,
            mode= 'markers',
        error_y= dict(
            array= prop_std,
            type= 'data',
            #symmetric= True,
            visible=True
        ),
            name= pop
        ))
        
        seg_counts= {
            sizes[x]: props[x] for x in range(len(props))
        }
        
        goods[pop]= seg_counts
    
    layout= go.Layout()
    Figure= go.Figure(data= fig, layout= layout)
    if plot_fig:
        iplot(Figure)
    
    return info_array,counts, goods



def varfac_collect(info_array,idx,counts,goods= {},func_plot= np.nanmean):
    '''
    factor function with optional plotly output and dictionary output of INNER FUNCTION, see Nseg_dict.
    '''
    #
    pop_dict= info_array[:,(idx-1)]
    pop_dict= dict_from_array(pop_dict)

    array_dict= {z: counts[g,:] for z,g in pop_dict.items()}
    
    fig= []
    
    
    for pop,pidx in pop_dict.items():
        #
        size_dict= info_array[pidx,idx]
        size_dict= np.array(size_dict,dtype= int)
        size_dict= dict_from_array(size_dict)
        
        sizes= sorted(size_dict.keys())[:-1]
        
        props= [func_plot(array_dict[pop][size_dict[z],:]) for z in sizes]

        seg_counts= {
            sizes[x]: props[x] for x in range(len(props))
        }

        goods[pop]= seg_counts
    
    
    return info_array,counts, goods




def size_la_plotly(info_array,idx,counts,goods= {},exc_max= True,func_plot= np.nanmean):
    '''
    factor local function w/ plotly output.
    '''
    #
    pop_dict= info_array[:,(idx-1)]
    pop_dict= dict_from_array(pop_dict)

    array_dict= {z: counts[g,:] for z,g in pop_dict.items()}
    
    fig= []

    for pop,pidx in pop_dict.items():
        #
        size_dict= info_array[pidx,idx]
        size_dict= np.array(size_dict,dtype= int)
        size_dict= dict_from_array(size_dict)
        
        sizes= sorted(size_dict.keys())[:-1]
        
        props= [func_plot(array_dict[pop][size_dict[z],:]) for z in sizes]
        fig.append(go.Scatter(
            x= sizes,
            y= props,
            mode= 'markers',
            name= pop
        ))
    
    layout= go.Layout()
    Figure= go.Figure(data= fig, layout= layout)
    iplot(Figure)
    
    return info_array,counts, goods

