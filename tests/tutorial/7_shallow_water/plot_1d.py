from scipy.optimize import fsolve
import math
import numpy as np
import h5py
import sys

import matplotlib
matplotlib.use('Pdf')
import matplotlib.pyplot as plt

import plot_settings as ps
import analytic_eq as eq

# plotting analytic solutions for height and velocity 
def analytic_fig(ax, time_l = [0,1,2,3], nx=320):
    x_range = np.linspace(-8,8, nx)
    oznacz = ['k', 'b', 'c', 'y', 'g', 'm', 'r']
    y0 = eq.d1_initial(x_range)
    for it, time in enumerate(time_l):
        lamb = eq.d1_lambda_evol(time)
        h = eq.d1_height(lamb, x_range)
        v = eq.d1_velocity(lamb, x_range)
        ax.plot(x_range, h, oznacz[it])
        ax.plot(x_range, v, oznacz[it]+ "--")
    ps.ticks_changes(ax)

#reading model output from text file and converting to an array
def reading_modeloutput(dirname, nt):
    f = h5py.File(dirname + "/timestep0000000" + str(nt) + ".h5", "r")
    h = np.array(f["h"])
    qx = np.array(f["qx"])
    f = h5py.File(dirname + "/coord.h5", "r")
    dt = f["T"].attrs["dt"][0]
    return h, qx, dt

#plotting together analytic solution and model output 
def analytic_model_fig(ax, x_range, h_m, v_m, time=1):
    lamb = eq.d1_lambda_evol(time)
    h_a = eq.d1_height(lamb, x_range)
    v_a = eq.d1_velocity(lamb, x_range)

    ax.plot(x_range, eq.d1_initial(x_range), 'k', x_range, h_a, 'b',
            x_range, h_m, "r")
    ax.plot(x_range, 0*x_range, "k-", x_range, v_a, 'b--',
            x_range, v_m, "r--")

    #ax.set_ylim(-2,2)
    ps.ticks_changes(ax)

# time_l - list of time levels for analytic solutions
# time - model time level used for comparison with analytic solution
# x_shift - shift between initial cond. in model and for analytic solution
def main(dir, casename_l, x_shift=8, time_l=[0,3], nt=300):
    plt.figure(1, figsize = (6,3))
    ax = plt.subplot(1,1,1)
    #plotting analytic solution
    analytic_fig(ax, time_l)
    plt.savefig("1d_analytic.pdf")
    #plotting comparison between analytic solution and model results for various options
    for ic, casename in enumerate(casename_l):
        plt.figure(ic+1, figsize = (6,3))
        ax = plt.subplot(1,1,1)

        print "plotting for " + casename + ", t/dt = " + str(nt)
        h_m, px_m, dt = reading_modeloutput(dir + casename, nt)
        # calculate velocity (momentum/height) only in the droplet region.
        v_m = np.where(h_m > 0,  px_m/h_m, 0)
                
        # choosing a plane of a cross section TODO: should be 160?
        # TODO: x_range/y_range should be calculated from hdf file!!
        analytic_model_fig(ax, np.linspace(-8,8,h_m.shape[0]), h_m, v_m, nt * dt)
 
        #ax.annotate(str(casename), xy=(0.01, 0.97), xycoords='axes fraction',
        #            fontsize=12, horizontalalignment='left', verticalalignment='top')
        plt.savefig(str(casename)+".pdf")

main("./", sys.argv[1:])
