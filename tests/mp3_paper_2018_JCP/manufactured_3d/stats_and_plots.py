import h5py
import numpy as np
import os, sys
sys.path.insert(1, os.path.join(os.path.dirname(__file__), '../python_scripts'))
from helpers import prepare_data, calc_convergence, opt2lab, nested_dict, save_conv
from conv_plot import conv_plot

def solution(geo_data, field_data, field, ny):
    time = 1.0
    nx, ny, nz = geo_data[ny]['g'].shape
    X = geo_data[ny]['di'] * np.arange(nx)
    Y = geo_data[ny]['dj'] * np.arange(ny)
    Z = geo_data[ny]['dk'] * np.arange(nz)
    X, Y, Z = np.meshgrid(X, Y, Z, indexing ='ij')
    solution = (2 + np.sin(X) * np.sin(time)) *\
               (2 + np.sin(Y) * np.sin(time)) *\
               (2 + np.sin(Z) * np.sin(time))
    return solution

def main():
    geo_data, field_data = prepare_data(sys.argv[1:])
    conv = calc_convergence(geo_data, field_data, '1.0', solution)

    plot_data = []
    norm = 'L2'
    field = 'psi'

    stat_file = open('manufactured_3d_conv.txt', 'w')
    save_conv(geo_data, conv, 'Manufactured solution in 3D', norm, field, stat_file)

    for opt in conv.keys():
        nys, errs = zip(*sorted(conv[opt][norm][field].items()))
        plot_data.append((nys, errs, opt))

    ord_data = []
    ord2 = lambda n : 3e-2 * (n / float(9)) ** (-2)
    ny2 = np.array([40, 140])
    nyt = 70
    ord_data.append((ny2, ord2(ny2), nyt, ord2(nyt+4), '2nd order', -92 - 180 / np.pi * np.arctan(-2)))

    ord3 = lambda n : 16e-3 * (n / float(9)) ** (-3)
    ny3 = np.array([40, 140])
    nyt = 70
    ord_data.append((ny3, ord3(ny3), nyt, ord3(nyt+4), '3rd order', -110 - 180 / np.pi * np.arctan(-3)))

    conv_plot(plot_data, ord_data, fname = 'manufactured_3d_conv.pdf')

main()
