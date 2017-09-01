import numpy as np
from helpers import set_rc_params
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FormatStrFormatter
import os
from analytical_formulae import asolution, pi

xlevels = np.linspace(0, 2 * pi, 17)
xlevels[-1] = 2 * pi - 0.01 
ylevels = np.linspace(-pi/2, pi / 2, 17)

def transform_iso_field(psi, X, Y, xo = 3 * pi / 2, yo = 0.):
    C = np.sin(yo) * np.sin(Y) + np.cos(yo) * np.cos(Y) * np.cos(X - xo)
    psi_i = np.where(C < 0, psi * np.nan, psi)
    return psi_i

def transform_iso_coord(X, Y, xo = 3 * pi / 2, yo = 0.):
    C = np.sin(yo) * np.sin(Y) + np.cos(yo) * np.cos(Y) * np.cos(X - xo)
    XI = np.cos(Y) * np.sin(X - xo)
    YI = np.cos(yo) * np.sin(Y) - np.sin(yo) * np.cos(Y) * np.cos(X - xo)

    XN = np.where(C < 0, X * np.nan, X)
    YN = np.where(C < 0, Y * np.nan, Y)
    XN = np.where(Y >= ylevels[-1], np.nan * XN, XN)
    return XI, YI, XN, YN

def panel_plot(geo_data, field_data, opt, time, ny, fname):
    set_rc_params(onlytex = True)

    nx, ny = geo_data[ny]['g'].shape
    X = geo_data[ny]['di'] * np.arange(nx)
    Y = geo_data[ny]['dj'] * (np.arange(ny) + 0.5) - pi / 2
    # extend Y values to cover the poles
    Y = np.concatenate([[-pi / 2], Y, [pi / 2]])
    X, Y = np.meshgrid(X, Y, indexing ='ij')

    init = asolution(0.0, X, Y)
    asol = asolution(float(time), X, Y)
    nsol = field_data[ny][opt][time]['psi']
    nsol = np.concatenate([np.zeros((nx, 1)), nsol, np.zeros((nx, 1))], axis = 1)
    errf = nsol - asol

    panel_fields = {'init' : init, 'asol' : asol, 'nsol' : nsol, 'errf' : errf}

    XI, YI, XN, YN = transform_iso_coord(X, Y)

    for pf in panel_fields:
        panel_fields[pf] = transform_iso_field(panel_fields[pf], X, Y)

    levels = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
    colors = ['white',
              'slateblue',
              'indigo',
              'mediumblue',
              'blue',
              'darkolivegreen',
              'olive',
              'lawngreen',
              'yellowgreen',
              'darkorange',
              'red',
              'firebrick'
              ]
    err_levels = [-0.021, -0.015, -0.009, -0.003,  0.003,  0.009, 0.015, 0.021]
    err_colors = [
              'slateblue',
              'blue',
              'indigo',
              'white',
              'darkorange',
              'red',
              'firebrick'
              ]
    fig, axarr = plt.subplots(2, 2, figsize = (15, 12))

    for ax in axarr.flatten():
        ax.axis('off')
        ax.contour(XI, YI, XN, levels = xlevels, linewidths = (0.5,), colors = ('k',), linestyles = ('dashed',))
        ax.contour(XI, YI, YN, levels = ylevels, linewidths = (0.5,), colors = ('k',), linestyles = ('dashed',))

    p = axarr[0, 0].contourf(XI, YI, panel_fields['init'], levels = levels, colors = colors)
    cb1 = fig.colorbar(p, ax = axarr[0, 0], ticks = levels[1:-1])

    p = axarr[0, 1].contourf(XI, YI, panel_fields['asol'], levels = levels, colors = colors)
    cb2 = fig.colorbar(p, ax = axarr[0, 1], ticks = levels[1:-1])

    p = axarr[1, 1].contourf(XI, YI, panel_fields['nsol'], levels = levels, colors = colors)
    cb3 = fig.colorbar(p, ax = axarr[1, 1], ticks = levels[1:-1])

    p = axarr[1, 0].contourf(XI, YI, panel_fields['errf'], levels = err_levels, colors = err_colors)
    cb4 = fig.colorbar(p, ax = axarr[1, 0], ticks = err_levels[1:-1])

    for cb in [cb1, cb2, cb3, cb4]:
        ticklabs = cb.ax.get_yticklabels()
        cb.ax.tick_params(labelsize=16) 
    
    axarr[0, 0].set_title('(a) Initial condition', y = 1.02, fontsize = 30)
    axarr[0, 1].set_title(r'(b) Exact solution $t = 12$ days', y = 1.02, fontsize = 30)
    axarr[1, 0].set_title(r'(c) Difference $t = 12$ days', y = 1.02, fontsize = 30)
    axarr[1, 1].set_title(r'(d) Numerical solution $t = 12$ days', y = 1.02, fontsize = 30)


    cb4.ax.set_yticklabels(ticklabs,ha='right')
    cb4.ax.yaxis.set_tick_params(pad=50) 
    plt.tight_layout(pad=1., w_pad=5., h_pad=1)
    plt.savefig(fname, transparent = True)
