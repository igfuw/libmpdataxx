import h5py
import numpy as np
import os, sys
sys.path.insert(1, os.path.join(os.path.dirname(__file__), '../mp3_paper_JCP/python_scripts'))
sys.path.insert(1, os.path.join(os.path.dirname(__file__), '../mp3_paper_JCP/reversing_deform'))
from mixing import calc_mixing_diags
from helpers import prepare_data, calc_convergence
from filaments import calc_filament_diags

def solution(geo_data, field_data, field, ny):
    key = list(field_data[ny].keys())[0]
    return field_data[ny][key]['0.0'][field]

def write_stats(conv, mixing_diags, filament_diags, opt, ny):
    outfile = open('stats_' + '_'.join(opt.decode('ASCII').split('|')) + '_i2_' + str(ny) + '.txt', 'w')

    errors_header = '{:6} {:8} {:8} {:8}\n'.format('field', 'L1', 'L2', 'Li')
    outfile.write(errors_header)

    for field in ['gh', 'cb', 'ccb', 'sc']:
        errs = {}
        for norm in ['L1', 'L2', 'Li']:
            errs[norm] = conv[opt][norm][field][ny]
        outfile.write('{:6} {L1:8.2e} {L2:8.2e} {Li:8.2e}\n'.format(field, **errs))

    filaments_header = '\n{:4}  {:6}\n'.format('tau', 'lf')
    outfile.write(filaments_header)
    for tau, lf in filament_diags[ny][opt]:
        outfile.write('{:4.2f}  {:5.1f}\n'.format(tau, lf))

    outfile.write('\nlr {:6.2e}\nlu {:6.2e}\nlo {:6.2e}\n'.format(*mixing_diags[ny][opt]))

def main():
    geo_data, field_data = prepare_data(sys.argv[1:])
    conv = calc_convergence(geo_data, field_data, '5.0', solution)

    nys = [120, 240]

    mixing_diags = calc_mixing_diags(geo_data, field_data, nys)
    filament_diags = calc_filament_diags(geo_data, field_data, nys)

    for ny in nys:
        for opt in conv.keys():
            write_stats(conv, mixing_diags, filament_diags, opt, ny)
main()
