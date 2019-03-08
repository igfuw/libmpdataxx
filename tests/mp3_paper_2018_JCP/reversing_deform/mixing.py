import numpy as np
from helpers import nested_dict, opt2lab, lab2pos, set_rc_params
import matplotlib
import matplotlib.pyplot as plt

def corr_func(x):
    return -0.8 * x ** 2 + 0.9

xmin = 0.1
xmax = 1.0
ymin = corr_func(xmin)
ymax = corr_func(xmax)

def line_func(x):
    a = (ymax - ymin) / (xmax - xmin)
    b = ymin - xmin * a
    return a * x + b

def root(x, y):
    sqrt_arg = 12 * (125 * y - 52) ** 3 + 29648025 * x ** 2
    if sqrt_arg < 0:
        print("problem, sqrt_arg < 0")
    c = 1. / 60 * (65340 * x + 12 * sqrt_arg ** 0.5) ** (1./3)
    root_uc = c + 1. / c * (13. / 75 - 5. / 12 * y)
    return min(max(xmin, root_uc), xmax)

def distance(x, y, root):
    xr = root
    yr = corr_func(root)
    return (((x - xr) / (xmax - xmin)) ** 2 + ((y - yr) / (ymax - ymin)) ** 2) ** 0.5

def save_mixing_diags(mixing_diags, ny, stat_file):
    ordered_opts = list(range(len(mixing_diags[ny].keys())))
    for opt in mixing_diags[ny].keys():
        ordered_opts[lab2pos[opt2lab[opt]]] = opt

    stat_file.write('Reversing deformational flow: mixing diagnostics for the N = {} ({} degree interval) simulations\n\n'
                    .format(ny, 180. / ny))
    stat_file.write('  ')
    for opt in ordered_opts:
        stat_file.write('{:>10}'.format(opt2lab[opt]))
    stat_file.write('\n')
    md_names = ['lr', 'lu', 'lo']
    for i in range(len(md_names)):
        stat_file.write('{:>2}'.format(md_names[i]))
        for opt in ordered_opts:
            stat_file.write('{:10.2e}'.format(mixing_diags[ny][opt][i]))
        stat_file.write('\n')

def calc_mixing_diags(geo_data, field_data, nys):
    mixing_diags = nested_dict()

    for ny in nys:
        g = geo_data[ny]['g']
        total_area = np.sum(g)

        for opt in field_data[ny]:
            cb, ccb = field_data[ny][opt]['2.5']['cb'].flatten(), field_data[ny][opt]['2.5']['ccb'].flatten()
            roots = map(root, cb , ccb)
            distances = map(distance, cb, ccb, roots)
            gdist = g.flatten() * list(distances)

            eps = 1e-7
            in_convex_hull = (ccb < corr_func(cb) + eps) & (ccb > line_func(cb) - eps)
            in_rectangle = (cb > xmin - eps) & (cb < xmax + eps) & (ccb > ymax - eps) & (ccb < ymin + eps)

            lr = np.sum(gdist[in_convex_hull]) / total_area
            lu = np.sum(gdist[in_rectangle & ~in_convex_hull]) / total_area
            lo = np.sum(gdist[~in_rectangle]) / total_area

            mixing_diags[ny][opt] = (lr, lu, lo)
    return mixing_diags

def plot_mixing(field_data, mixing_diags, ny, fname = 'mixing.pdf'):
    set_rc_params()
    matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

    plt.clf()
    fig, axarr = plt.subplots(2, 2, sharex='col', sharey = 'row', figsize=(14, 12))

    c = 'red'
    s = 1
    marker = '.'
    edgecolor = 'r'
    fs = 22
    lw = 2

    x = np.linspace(xmin, xmax)

    for opt in field_data[ny]:
        cb, ccb = field_data[ny][opt]['2.5']['cb'].flatten(), field_data[ny][opt]['2.5']['ccb'].flatten()

        lab = opt2lab[opt]

        lab2num = {'Mp3' : (0, 1), 'Mp3cc' : (0, 0), 'Mg2No' : (1, 0), 'Mg3No' : (1, 1)}
        num = lab2num[lab]

        axarr[num].scatter(cb, ccb, c = c, s = s, marker = marker, edgecolor = edgecolor, rasterized = True)
        axarr[num].plot(x, corr_func(x), 'k-', lw = lw)
        axarr[num].plot(x, ymin * x / x, 'k-', lw = lw)
        axarr[num].plot(xmax * x / x, np.linspace(ymax, ymin), 'k-', lw = lw)
        axarr[num].plot(x, line_func(x), 'k-', lw = lw)

        axarr[num].set_ylim([0, 1])
        axarr[num].set_xlim([0, 1.1])
        axarr[num].set_title(lab, fontsize = fs + 8, y = 1.02)
        if (num not in [(0, 0), (0, 1)]):
            axarr[num].set_xlabel('$\\chi$', fontsize = fs+2)
        if (num not in [(1, 1), (0, 1)]):
            axarr[num].set_ylabel('$\\xi$', fontsize = fs+2)

        axarr[num].xaxis.set_tick_params(direction = 'in', top = True, labelsize = 20)
        axarr[num].xaxis.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        #if (num in [(1, 0), (1, 1)]):
        axarr[num].yaxis.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        axarr[num].yaxis.set_tick_params(direction = 'in', right = True, labelsize = 20)
        axarr[num].minorticks_off()

        axarr[num].grid()

        lr, lu, lo = mixing_diags[ny][opt]

        expformat = '{} \\times 10^{{{}}}'

        lrs = '{:.2e}'.format(lr)
        l_r_d, l_r_e = (lrs[0:4], int(lrs[5:]))

        lus = '{:.2e}'.format(lu)
        l_u_d, l_u_e = lus[0:4], int(lus[5:])

        if lo != 0:
            los = '{:.2e}'.format(lo)
            lo = (los[0:4], int(los[5:]))
            loformat = expformat
        else:
            lo = (lo,)
            loformat = '{:g}'

        test2 = (
                '''\\begin{{{{align}}}}'''
                '''l_r &= {l_r_format}\\\\'''
                '''l_u &= {l_u_format}\\\\'''
                '''l_o &= {l_o_format}\\\\'''
                '''\\end{{{{align}}}}'''
               )
        test = test2.format(l_r_format = expformat,
                            l_u_format = expformat,
                            l_o_format = loformat)

        axarr[num].text(0.1,
                        0.05,
                        test.format(l_r_d, l_r_e, l_u_d, l_u_e, *lo),
                        size = 28,
                        multialignment = 'left'
                       )
    plt.tight_layout(pad=1., w_pad=5., h_pad=1)
    plt.savefig(fname)
