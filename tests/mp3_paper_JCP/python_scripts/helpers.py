import h5py
import numpy as np
from sys import argv
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict

# https://stackoverflow.com/questions/5369723/multi-level-defaultdict-with-variable-depth/8702435#8702435
nested_dict = lambda: defaultdict(nested_dict)

opt2lab = {
            'nug|abs|div_2nd|div_3rd' : 'Mp3',
            'nug|div_2nd|div_3rd' : 'Mp3',
            'nug|abs|dfl|tot' : 'Mp3cc',
            'nug|dfl|tot' : 'Mp3cc',
            'nug|abs|tot' : 'Mp3cc',
            'nug|tot' : 'Mp3cc',
            'nug|iga|div_2nd|div_3rd|fct' : 'Mg3No',
            'nug|iga|fct' : 'Mg2No',
            'nug|iga|dfl|fct' : 'Mg2No'
          }

lab2pos = {'Mp3' : 0, 'Mp3cc' : 1, 'Mg3No' : 3, 'Mg2No' : 2}

def prepare_data(dirnames, chosen_times = []):
    geo_data = nested_dict()
    field_data = nested_dict()

    for dirname in dirnames:
        fnames = sorted(filter(lambda s : s.endswith('.h5'), os.listdir(dirname)))
        fnames = [dirname + '/' + name for name in fnames]

        const_file = h5py.File(fnames[0], 'r')

        g = const_file['G'][:,:]
        ndims = len(g.shape) 
        ny = g.shape[1]

        geo_data[ny] = {'g' : g}
        geo_data[ny]['di'] = const_file['advection'].attrs['di'][0]
        geo_data[ny]['dj'] = const_file['advection'].attrs['dj'][0]
        if ndims > 2:
            geo_data[ny]['dk'] = const_file['advection'].attrs['dk'][0]
        
        opts = const_file['advection'].attrs['opts'][0].decode('ASCII')
        times = const_file['T'][:]
      
        if len(chosen_times) > 0:
            indices_t = []
            for m in range(len(times)):
                for ct in chosen_times:
                    if np.isclose(times[m], float(ct)):
                        indices_t.append(m)
        else:
            indices_t = range(len(times))

        for t in indices_t:
            ts = '{:.1f}'.format(times[t])
            data_file = h5py.File(fnames[t+1], 'r')
            for field in data_file['/']:
                field_data[ny][opts][ts][field] = data_file[field][:,:]
        
    return geo_data, field_data

def calc_convergence(geo_data, field_data, at_time, solution):
    norms = [
             ('L1', lambda s, n, g : np.sum(g * abs(n - s)) / np.sum(g * np.abs(s))),
             ('L2', lambda s, n, g : np.sqrt(np.sum(g * (n - s) ** 2) / np.sum(g * s ** 2))),
             ('Li', lambda s, n, g : np.max(np.abs(n - s)) / np.max(np.abs(s)))
            ]
    conv = nested_dict()
    for ny in field_data.keys():
        g = geo_data[ny]['g']
        for opt in field_data[ny].keys():
            # 0.0 time always available
            for field in field_data[ny][opt]['0.0'].keys():
                sol = solution(geo_data, field_data, field, ny)
                for n_name, n_func in norms:
                    conv[opt][n_name][field][ny] = n_func(sol, field_data[ny][opt][at_time][field], g)
    return conv

def set_rc_params(tickfs = 35, onlytex = False):
    mpl.rcdefaults()
    mpl.rc('text', usetex=True)
    mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    if not onlytex:
        mpl.rcParams['xtick.labelsize'] = tickfs
        mpl.rcParams['ytick.labelsize'] = tickfs
        mpl.rcParams['legend.fontsize'] = 45
        tw = 3
        ts = 12
        plt.rc('xtick.major', width = tw)
        plt.rc('xtick.major', size = ts)
        plt.rc('ytick.major', width = tw)
        plt.rc('ytick.major', size = ts)
        plt.rc('axes', linewidth = 3)

def save_conv(geo_data, conv, test, norm, field, stat_file):
    stat_file.write('{}: convergence in the {} error norm\n\n'.format(test, norm))
    stat_file.write('{:>3}'.format('N'))

    ordered_opts = list(range(len(conv.keys())))
    for opt in conv.keys():
        ordered_opts[lab2pos[opt2lab[opt]]] = opt

    for opt in ordered_opts:
        stat_file.write('{:>10}'.format(opt2lab[opt]))
    stat_file.write('\n')

    for ny in sorted(geo_data.keys()):
        stat_file.write('{:3d}'.format(ny))
        for opt in ordered_opts:
            stat_file.write('{:10.2e}'.format(conv[opt][norm][field][ny]))
        stat_file.write('\n')
