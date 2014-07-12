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

def d2_velocity(lamb, x, y):
    lamb_t = 2**0.5 * (1 - lamb**-2)**0.5
    return np.where(d2_rad2(x,y) <= lamb**2, x * lamb_t / lamb, 0)

#initial condition
def d2_initial(x,y):
    return np.where(d2_rad2(x,y)<=1, 1-d2_rad2(x,y), 0)

# for elliptic:
def d2_el_height(lamb_x, lamb_y, x, y):
    return np.where(x**2/lamb_x**2 + y**2/lamb_y**2 <= 1.,
                    1./lamb_x/lamb_y * (1. - x**2/lamb_x**2 - y**2/lamb_y**2), 0.)

#can be used for both - x and y
def d2_el_velocity(lamb, lamb_t, lamb_other, x, x_other):
    return np.where(x**2/lamb**2 + x_other**2/lamb_other**2 <= 1.,
                    x * lamb_t / lamb, 0)

#for derivation TODO! - rewrite
def deriv(y,t):
    # return derivatives of [lambda_x, dlambda_x/dt, lambda_y, dlambda_y/dt
    return np.array([y[1], 2. / y[0]**2 / y[2], y[3], 2. / y[0] / y[2]**2])


def d2_el_lamb_lamb_t_evol(time, lamb_x0, lamb_y0):
    time_lin = np.linspace(0.0,time,101)
    yinit = np.array([lamb_x0, 0., lamb_y0, 0.]) # initial values (velocity is 0.)
    print "lamb_x0, lamb_y0", lamb_x0, lamb_y0, time
    y = odeint(deriv,yinit,time_lin)
    print "time, lambda_x, dlambda_x/dt, lambda_y, dlambda_y/dt", time_lin[-1], y[-1,:]
    return y[-1,:]
