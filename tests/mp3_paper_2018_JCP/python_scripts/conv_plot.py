import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FormatStrFormatter
import os
from helpers import set_rc_params, opt2lab, lab2pos

def conv_plot(plot_data, ord_data, fname = 'fig.pdf'):
    set_rc_params()
    fs = 45

    lab2ls = {'Mp3' : 'ks-', 'Mp3cc' : 'bo-', 'Mg3No' : 'gD-', 'Mg2No' : 'r^-'}
    lab2dsh = {'Mp3' : [], 'Mp3cc' : [6, 1], 'Mg3No' : [1, 1], 'Mg2No' : [6, 1, 1 ,1]}

    fig, axarr = plt.subplots(1, 1, figsize=(14,12))

    plots = {}
    for x, y, opt in plot_data:
        zo = 2.0
        lab = opt2lab[opt]
        if lab == 'Mg3No':
            zo = 2.1
        p, = axarr.loglog(x, y, lab2ls[lab], label = lab, lw = 6, ms = 16, dashes = lab2dsh[lab], zorder = zo)
        plots[lab] = p

    for x, y, x0, y0, t, a in ord_data:
        axarr.loglog(x, y, color = 'gray', lw = 6)
        axarr.text(x0, y0, t, rotation = a, size = 28, backgroundcolor = 'w', zorder=2.1)

    axarr.xaxis.set_ticks(plot_data[0][0])
    axarr.xaxis.set_tick_params(direction = 'in')
    axarr.yaxis.set_tick_params(direction = 'in')
    axarr.minorticks_off()

    axarr.get_xaxis().set_major_formatter(ScalarFormatter())

    axarr.set_xlabel('$N$', fontsize = fs)
    axarr.set_ylabel('$\ell_2$ error norm', fontsize = fs)

    labels = list(range(len(plots.keys())))
    for lab in plots.keys():
        labels[lab2pos[lab]] = lab
    handles = [plots[lab] for lab in labels]
    l = axarr.legend(handles, labels, fontsize=fs, numpoints=1, edgecolor = 'k', framealpha = 1.0, loc = 3, handlelength = 2.5)
    l.get_frame().set_linewidth(3)

    axarr.grid(color = 'k', lw = 2, linestyle = ':')

    plt.tight_layout(pad=1., w_pad=5., h_pad=1)
    plt.savefig(fname)
