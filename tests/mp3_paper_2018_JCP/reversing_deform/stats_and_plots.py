import h5py
import numpy as np
import os, sys
sys.path.insert(1, os.path.join(os.path.dirname(__file__), '../python_scripts'))
from helpers import prepare_data, calc_convergence, save_conv
from conv_plot import conv_plot
from panel_plot import panel_plot
from time import clock
from mixing import calc_mixing_diags, save_mixing_diags, plot_mixing

def solution(geo_data, field_data, field, ny):
    key = list(field_data[ny].keys())[0]
    return field_data[ny][key]['0.0'][field]

def main():
    geo_data, field_data = prepare_data(sys.argv[1:], chosen_times = ['0.0', '2.5', '5.0'])
    conv = calc_convergence(geo_data, field_data, '5.0', solution)

    plot_data = []
    norm = 'L2'
    field = 'gh'

    for opt in conv.keys():
        nys, errs = zip(*sorted(conv[opt][norm]['gh'].items()))
        plot_data.append((nys, errs, opt))

    stat_file = open('reversing_deform_conv.txt', 'w')
    save_conv(geo_data, conv, 'Reversing deformational flow', norm, field, stat_file)
    stat_file.close()

    ord_data = []
    ord2 = lambda n : 1.5e-1 * (n / 120.) ** (-2)
    ny2 = np.array([400, 1000])
    nyt = 580
    ord_data.append((ny2, ord2(ny2), nyt, ord2(nyt+10), '2nd order', -95 - 180 / np.pi * np.arctan(-2)))

    ord3 = lambda n : 1e-1 * (n / 120.) ** (-3)
    ny3 = np.array([400, 1000])
    nyt = 640
    ord_data.append((ny3, ord3(ny3), nyt, ord3(nyt+20), '3rd order', -115 - 180 / np.pi * np.arctan(-3)))

    conv_plot(plot_data, ord_data, fname = 'reversing_deform_conv.pdf')

    panel_plot(field_data[120]['nug|iga|div_2nd|div_3rd|fct']['2.5'], fname = 'reversing_deform_panel.pdf')
   
    mixing_ny = 240
    mixing_diags = calc_mixing_diags(geo_data, field_data, nys = [mixing_ny])

    stat_file = open('reversing_deform_mixing.txt', 'w')
    save_mixing_diags(mixing_diags, mixing_ny, stat_file)
    stat_file.close()

    plot_mixing(field_data, mixing_diags, mixing_ny, fname = 'reversing_deform_mixing.pdf')

main()
