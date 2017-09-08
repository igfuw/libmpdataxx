import numpy as np

pi = np.pi

def xpf(x, y, x0, y0):
    ret = np.arctan2(np.cos(y) * np.sin(x - x0), np.cos(y) * np.sin(y0) * np.cos(x - x0) - np.cos(y0) * np.sin(y));
    return np.where(ret <= 0, ret + 2 * pi, ret)

def ypf(x, y, x0, y0):
    return np.arcsin(np.sin(y) * np.sin(y0) + np.cos(y) * np.cos(y0) * np.cos(x - x0))

def ixpf(x, y, x0, y0):
    ret = x0 + np.arctan2(np.cos(y) * np.sin(x), np.sin(y) * np.cos(y0) + np.cos(y) * np.cos(x) * np.sin(y0));
    return np.where(ret <= 0, ret + 2 * pi, ret)

def iypf(x, y, x0, y0):
    return np.arcsin(np.sin(y) * np.sin(y0) - np.cos(y) * np.cos(y0) * np.cos(x));

def asolution(t, x, y, a = pi / 2):
    x0 = 3 * pi / 2
    y0 = 0.0
    u0 = 2 * pi / 12
    v0 = 2 * pi / 12

    xtmp = xpf(x0, y0, pi, pi / 2 - a)
    xtmp += u0 * t
    ytmp = ypf(x0, y0, pi, pi / 2 - a)

    xc = ixpf(xtmp, ytmp, pi, pi / 2 - a)
    yc = iypf(xtmp, ytmp, pi, pi / 2 - a)

    r = 3 * np.cos(ypf(x, y, xc, yc))
    omg = np.where(r != 0, u0 * 3 * np.sqrt(2.) / (2 * r) * np.tanh(r) / (np.cosh(r) ** 2), 0)
    solution = 1 - np.tanh(r / 5 * np.sin(xpf(x, y, xc, yc) - omg * t))
    return solution
