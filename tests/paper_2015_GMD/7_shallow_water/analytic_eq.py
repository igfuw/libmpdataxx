import numpy as np
import math
from scipy.optimize import fsolve
from scipy.integrate import odeint

#eq. 5.6 (Schar&Smolarkiewicz, 1996)
def d1_lambda_eq(x, *time):
    return 0.5 * ((x*(x-1.))**0.5 + math.log((x-1.)**0.5 + x**0.5)) - time

# finging roots of lambda_eq
def d1_lambda_evol(time, x0=1):
    return fsolve(d1_lambda_eq, x0, args=time)

#eq. 5.4 (Schar&Smolarkiewicz, 1996)
def d1_height(lamb, x):
    return np.where(x**2 <= lamb**2, lamb**-1 * (1. - (x/lamb)**2), 0)

#eq. 5.5 (Schar&Smolarkiewicz, 1996)
def d1_velocity(lamb, x):
    lamb_t = 2 * (1 - lamb**-1)**0.5
    return np.where(x**2 <= lamb**2, x * lamb_t / lamb, 0)

#eq. 5.3 (Schar&Smolarkiewicz, 1996)
def d1_initial(x):
    return np.where(x**2<=1, 1-x**2, 0)


def d2_rad2(x,y):
    return x**2 + y**2

def d2_lambda_evol(time):
    return (2*time**2 + 1)**0.5

def d2_height(lamb, x, y):
    return np.where(d2_rad2(x,y) <= lamb**2, lamb**-2 * (1. - d2_rad2(x,y)/lamb**2), 0)

# height in 2d
def d2_height_plane(lamb, x, y):
    X, Y = np.meshgrid(x, y)
    return np.where(d2_rad2(X,Y) <= lamb**2, lamb**-2 * (1. - d2_rad2(X,Y)/lamb**2), 0)

def d2_velocity(lamb, x, y):
    lamb_t = 2**0.5 * (1 - lamb**-2)**0.5
    return np.where(d2_rad2(x,y) <= lamb**2, x * lamb_t / lamb, 0)

#initial condition
def d2_initial(x,y):
    return np.where(d2_rad2(x,y)<=1, 1-d2_rad2(x,y), 0)
