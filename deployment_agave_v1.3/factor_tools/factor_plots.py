
import os
import numpy as np

import matplotlib.pyplot as plt
import plotly
#import chart_studio.plotly as py
import plotly.graph_objs as go
from plotly.graph_objs import *
from plotly.offline import iplot

## plot functions

def plot_power_3D(rate_stats,rate_range,func_select= np.nanmean,Nsteps_tuple= (4,4),
                  xtitle= '', ytitle= '', ztitle= ''):
    """
    """
    bi_select= list(rate_stats.keys(),)
    
    fig= []
    d= 0
    
    colorscale= ['rgb(100,0,0)','rgb(100,0,100)','rgb(0,100,0)']
    colorscale= [[x / len(colorscale),colorscale[x]] for x in range(len(colorscale))]
    print(colorscale)
    for rate in rate_range:

        props= [(*bi,func_select(rate_stats[bi][rate]['pval'])) for bi in bi_select]
        
        props= np.array(props)
        #
        color_grad= np.zeros(Nsteps_tuple)
        color_grad+= d / len(rate_range)
        #
        fig.append(go.Surface(
            x= props[:,0].reshape(Nsteps_tuple),
            y= props[:,1].reshape(Nsteps_tuple),
            z= props[:,2].reshape(Nsteps_tuple),
            cmin=0,
            cmax=1,
            cauto= False,
            showscale=False, opacity=1,
            surfacecolor= color_grad,
            colorscale= colorscale,
            name= 'rate fold: {}'.format(rate)
        ))
        d += 1
    
    layout= go.Layout(
        showlegend=True,
        scene= go.layout.Scene(
            xaxis= go.layout.scene.XAxis(title= xtitle),
            yaxis= go.layout.scene.YAxis(title= ytitle),
            zaxis= go.layout.scene.ZAxis(title= ztitle)
        ),
        height= 800,
        width= 900
    )
    
    Figure= go.Figure(data= fig, layout= layout)
    #Figure.show()
    iplot(Figure)



def plot_suraceCt(rate_stats,rate,func_select= np.nanmean,
                 height= 800,width= 800,title_add= '', xtitle= 's1', ytitle= 's2',
                 savefig= ''):
    '''
    plot surface contour.
    '''
    bi_select= list(rate_stats.keys())
    
    fig= []
    d= 0

    props= [(*bi,func_select(rate_stats[bi][rate]['pval'])) for bi in bi_select]    
    props= np.array(props)
    gradient= props[:,2]
    
    fig.append(go.Contour(
        x= props[:,0],
        y= props[:,1],
        z= gradient
    ))
    d += 1
    
    layout= go.Layout(
        title= 'rate: {} '.format(rate) + title_add,
        height= height,
        width= width,
        xaxis= dict(
            title= xtitle
        ),
        yaxis= dict(
            title= ytitle
        )
    )
    
    Figure= go.Figure(data= fig, layout= layout)

    iplot(Figure)




import matplotlib.cm as cm

def plot_contourPLT(rate_stats,rate,func_select= np.nanmean,
                 height= 10,width= 15,title_add= '', xtitle= 's1', ytitle= 's2',
                 savefig= '', xlim= [], ylim= [], clim= [], levels= 6):
    '''
    plot surface contour.
    '''
    bi_select= list(rate_stats.keys())
    
    fig= []
    d= 0
    
    props= [(*bi,func_select(rate_stats[bi][rate]['pval'])) for bi in bi_select]    
    props= np.array(props)
    
    gradient= props[:,2]
    
    plt.inferno()
    fig=plt.figure(figsize=(width,height))

    X= props[:,0]
    Y= props[:,1]
    
    dims= [len(set(X)),len(set(Y))]
    dims= tuple(dims)
    
    X= X.reshape(dims)
    Y= Y.reshape(dims)
    gradient= gradient.reshape(dims)
    
    if clim:
        ctr_levels= np.linspace(clim[0],clim[1],levels+1)
        cp = plt.contourf(X, Y, gradient, levels= ctr_levels, cmap=cm.inferno)
    else:
        cp = plt.contourf(X, Y, gradient,levels= levels,cmap= cm.inferno)
    
    fig.colorbar(cp)
    plt.title('rate: {} '.format(rate) + title_add)
    plt.xlabel(xtitle)
    if xlim:
        plt.xlim(*xlim)
    plt.ylabel(ytitle)
    if ylim:
        plt.ylim(*ylim)
    if savefig:
        os.makedirs(os.path.dirname(savefig), exist_ok=True)
        plt.savefig(savefig)
        plt.close()
    
    else:
        plt.show()
        plt.close()



def plot_fixedSizes(comb_dict, height= 800, width= 800, title= ''):
    '''
    plot comb_dict
    '''

    samp_max= max([max(g['sizes']) for g in comb_dict.values()])
    fig= [go.Scatter(
        x= comb_dict[i]['sizes'],
        y= comb_dict[i]['rates'],
        mode= 'lines',
        name= '-'.join(i)
    ) for i in comb_dict.keys()]

    layout= go.Layout(
        title= title,
        xaxis= dict(
            range= [-1,samp_max+1],
            title= 'Nsamp - same'
        ),
        yaxis= dict(
          title= 'min rate fold'  
        ),
        height= height,
        width= width
    )

    Figure= go.Figure(data= fig, layout= layout)
    iplot(Figure)
