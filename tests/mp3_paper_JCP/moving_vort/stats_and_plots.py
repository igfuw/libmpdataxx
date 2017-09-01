import h5py
import numpy as np
import os, sys
sys.path.insert(1, os.path.join(os.path.dirname(__file__), '../python_scripts'))
from helpers import prepare_data, calc_convergence, save_conv
from conv_plot import conv_plot
from time import clock
from analytical_formulae import asolution, pi
from panel_plot import panel_plot

def solution(geo_data, field_data, field, ny):
    nx, ny = geo_data[ny]['g'].shape
    X = geo_data[ny]['di'] * np.arange(nx)
    Y = geo_data[ny]['dj'] * (np.arange(ny) + 0.5) - pi / 2
    X, Y = np.meshgrid(X, Y, indexing ='ij')
    
    return asolution(12.0, X, Y)

def main():
    geo_data, field_data = prepare_data(sys.argv[1:], chosen_times = ['0.0', '12.0'])
    conv = calc_convergence(geo_data, field_data, '12.0', solution)

    plot_data = []
    norm = 'L2'
    field = 'psi'

    stat_file = open('moving_vort_conv.txt', 'w')
    save_conv(geo_data, conv, 'Moving vortices', norm, field, stat_file)

    for opt in conv.keys():
        nys, errs = zip(*sorted(conv[opt][norm][field].items()))
        plot_data.append((nys, errs, opt))

    ord_data = []
    ord2 = lambda n : 2e-0 * (n / 24.) ** (-2)
    ny2 = np.array([300, 900])
    nyt = 450
    ord_data.append((ny2, ord2(ny2), nyt, ord2(nyt+10), '2nd order', -94 - 180 / pi * np.arctan(-2)))

    ord3 = lambda n : 1e-1 * (n / 24.) ** (-3)
    ny3 = np.array([300, 900])
    nyt = 480
    ord_data.append((ny3, ord3(ny3), nyt, ord3(nyt+20), '3rd order', -113 - 180 / pi * np.arctan(-3)))

    conv_plot(plot_data, ord_data, fname = 'moving_vort_conv.pdf')

    panel_plot(geo_data, field_data, opt = 'nug|abs|div_2nd|div_3rd', time = '12.0', ny = 192, fname = 'moving_vort_panel.pdf')

main()
