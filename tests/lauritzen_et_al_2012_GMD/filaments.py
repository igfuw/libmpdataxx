import numpy as np
import matplotlib
matplotlib.use('Pdf')
import matplotlib.pyplot as plt

def calc_filament_diags(g, field_data):
    taus = np.linspace(0.10, 1.0, 19)
    lfs = []
    eps = 1e-12
    for tau in taus:
        area_0 = np.sum(g[field_data['cb']['initial'] >= tau - eps])
        area_h = np.sum(g[field_data['cb']['halfway'] >= tau - eps])
        if area_0 < eps:
            lf = 0
        else:
            lf = 100 * area_h / area_0
        lfs.append(lf)
    return zip(taus, lfs)

def plot_filament_diags(dirname, filament_diags, deg):
    plt.clf()
    taus, lfs = zip(*filament_diags)
    plt.plot(taus, lfs, 'ks-', lw = 2, ms = 8)
    fs = 22
    plt.xlabel(r'$\tau$', fontsize = fs)
    plt.ylabel(r'$l_f(\tau, T / 2)$', fontsize = fs)
    plt.ylim([0, 140])
    plt.title('filament preservation diagnostics at ${}\\degree$'.format(deg), fontsize = fs)
    plt.grid()
    plt.savefig('filaments_' + dirname + '.pdf')
