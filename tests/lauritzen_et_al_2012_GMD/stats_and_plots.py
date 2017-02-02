import h5py
import numpy as np
from sys import argv
import os

from filaments import calc_filament_diags, plot_filament_diags
from mixing import calc_mixing_diags, plot_mixing

def prepare_data(dirname):
    fnames = sorted(filter(lambda s : s.endswith('.h5'), os.listdir(dirname)))
    fnames = [dirname + '/' + name for name in fnames]
    
    fields = ['gh', 'cb', 'sc', 'ccb']
    states = ['initial', 'halfway', 'final']

    field_data = {}
    g = h5py.File(fnames[0], 'r')['G'][:,:]
    for f in fields:
        field_data[f] = {}
        for i in range(len(states)):
            data_file = h5py.File(fnames[i+1], 'r')
            field_data[f][states[i]] = data_file[f][:,:]
    return (g, field_data)

def calc_errors(g, field_data):
    norms = [
             ('L1', lambda i, f : np.sum(g * abs(f - i)) / np.sum(g * np.abs(i))),
             ('L2', lambda i, f : np.sqrt(np.sum(g * (f - i) ** 2) / np.sum(g * i ** 2))),
             ('Li', lambda i, f : np.max(np.abs(f - i)) / np.max(np.abs(i)))
            ]
    errors = {}
    for f in field_data.keys():
        errors[f] = {}
        for n_name, n_func in norms:
            errors[f][n_name] = n_func(field_data[f]['initial'], field_data[f]['final'])
    return errors

def write_stats(dirname, errors, filament_diags, mixing_diags):
    outfile = open('stats_' + dirname + '.txt', 'w')

    errors_header = '{:6} {:9} {:9} {:9}\n'.format('field', 'L1', 'L2', 'Li')
    outfile.write(errors_header)
    for f in ['gh', 'cb', 'ccb', 'sc']:
        outfile.write('{:6} {L1:9.3e} {L2:9.3e} {Li:9.3e}\n'.format(f, **errors[f]))

    filaments_header = '\n{:4}  {:6}\n'.format('tau', 'lf')
    outfile.write(filaments_header)
    for tau, lf in filament_diags:
        outfile.write('{:4.2f}  {:5.1f}\n'.format(tau, lf))

    outfile.write('\nlr {:6.4e}\nlu {:6.4e}\nlo {:6.4e}'.format(*mixing_diags))

def main():
    for dirname in argv[1:]:
        g, field_data = prepare_data(dirname)
        errors = calc_errors(g, field_data)
        filament_diags = calc_filament_diags(g, field_data)
        mixing_diags = calc_mixing_diags(g, field_data)

        deg = str(180. / g.shape[1])
        plot_filament_diags(dirname, filament_diags, deg)
        plot_mixing(dirname, field_data, deg)
        write_stats(dirname, errors, filament_diags, mixing_diags)

main()
