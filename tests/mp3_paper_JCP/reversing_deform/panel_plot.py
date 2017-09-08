import numpy as np
from helpers import set_rc_params
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FormatStrFormatter
import os

def panel_plot(plot_data, fname = 'panel.pdf'):
    set_rc_params(tickfs = 30)

    levels = np.linspace(-0.05, 1.15, 25)

    fig, axarr = plt.subplots(2, 2, figsize=(20, 14), sharex = 'col', sharey = 'row')
    field2pos = {'gh' :  (0, 0), 'cb' : (0, 1), 'ccb' : (1, 1), 'sc' : (1, 0) }
    field2title = {'gh' :  '(a) Gaussian hills',
                   'cb' : '(b) cosine bells',
                   'ccb' : "(d) 'correlated' cosine bells",
                   'sc' : '(c) slotted cylinders'}

    for f in plot_data.keys():
        field = np.transpose(plot_data[f])
        ny, nx = field.shape
        cmap = plt.get_cmap('terrain')
        i = field2pos[f]
        p = axarr[i].contourf(field, levels = levels, zorder = 1, cmap = cmap)
        axarr[i].contour(field, levels = levels, zorder = 1, cmap = cmap)
        axarr[i].set_title(field2title[f], fontsize = 30, y = 1.02)
    
    for ax in axarr.flatten():
        ax.set_xticks([0, (nx - 1) / 4, (nx - 1) / 2, 3 * (nx - 1) / 4, nx - 1])
        ax.set_xticklabels([r'$0$', r'$\pi$/2', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
        
        ax.set_yticks([0, (ny+0.5) / 2, ny - 1])
        ax.set_yticklabels([r'$-\pi/2$', r'$0$', r'$\pi/2$'])
        ax.grid(color = 'w', lw = 2, linestyle = ':')

    plt.subplots_adjust(bottom=0.2)

    cbaxes = fig.add_axes([0.1, 0.1, 0.8, 0.02]) 
    plt.colorbar(p, cax = cbaxes, ticks = levels[1:-1:2], orientation = 'horizontal')

    #plt.tight_layout(pad=1., w_pad=5., h_pad=1)
    plt.savefig(fname, bbox_inches='tight')
