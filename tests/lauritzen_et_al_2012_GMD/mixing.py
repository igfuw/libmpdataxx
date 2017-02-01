import numpy as np
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
        print "problem, sqrt_arg < 0"
    c = 1. / 60 * (65340 * x + 12 * sqrt_arg ** 0.5) ** (1./3)
    root_uc = c + 1. / c * (13. / 75 - 5. / 12 * y)
    return min(max(xmin, root_uc), xmax)

def distance(x, y, root):
    xr = root
    yr = corr_func(root)
    return (((x - xr) / (xmax - xmin)) ** 2 + ((y - yr) / (ymax - ymin)) ** 2) ** 0.5

def calc_mixing_diags(g, field_data):
    cb, ccb = field_data['cb']['halfway'].flatten(), field_data['ccb']['halfway'].flatten()
    roots = map(root, cb , ccb)
    distances = map(distance, cb, ccb, roots)
    total_area = np.sum(g)
    gdist = g.flatten() * distances

    eps = 1e-7
    in_convex_hull = (ccb < corr_func(cb) + eps) & (ccb > line_func(cb) - eps)
    in_rectangle = (cb > xmin - eps) & (cb < xmax + eps) & (ccb > ymax - eps) & (ccb < ymin + eps)

    lr = np.sum(gdist[in_convex_hull]) / total_area
    lu = np.sum(gdist[in_rectangle & ~in_convex_hull]) / total_area
    lo = np.sum(gdist[~in_rectangle]) / total_area
    return lr, lu, lo

def plot_mixing(dirname, field_data, deg):
    plt.clf()
    cb, ccb = field_data['cb']['halfway'].flatten(), field_data['ccb']['halfway'].flatten()
    c = 'red' 
    s = 1
    marker = '.'
    edgecolor = 'r'
    fs = 22
    lw = 2

    plt.scatter(cb, ccb, c = c, s = s, marker = marker, edgecolor = edgecolor, rasterized = True)
    x = np.linspace(xmin, xmax)
    plt.plot(x, corr_func(x), 'k-', lw = lw)
    plt.plot(x, ymin * x / x, 'k-', lw = lw)
    plt.plot(xmax * x / x, np.linspace(ymax, ymin), 'k-', lw = lw)
    plt.plot(x, line_func(x), 'k-', lw = lw)

    plt.ylim([0, 1])
    plt.xlim([0, 1.1])
    plt.xlabel('cb', fontsize = fs)
    plt.ylabel('ccb', fontsize = fs)
    plt.title('mixing at ${}\\degree$'.format(deg), fontsize = fs)
    plt.grid()

    plt.savefig('mixing_' + dirname + '.pdf')
