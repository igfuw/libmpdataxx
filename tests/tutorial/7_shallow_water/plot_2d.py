from scipy.optimize import fsolve
import math
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib.ticker as ticker
import numpy as np
import h5py
import sys
import plot_settings as ps
import analytic_eq as eq

# plotting analytic solutions for height and velocity 
def analytic_fig(ax, time_l = [0,1,2,3], x_range = np.linspace(-8,8,320),
                              y_range = np.zeros(320)):
    oznacz = ['k', 'b', 'c', 'y', 'g', 'm', 'r']
    y0 = eq.d2_initial(x_range, y_range)

    for it, time in enumerate(time_l):
        lamb = eq.d2_lambda_evol(time)
        h = eq.d2_height(lamb, x_range, y_range)
        v = eq.d2_velocity(lamb, x_range, y_range)
        ax.plot(x_range, h, oznacz[it])
        ax.plot(x_range, v, oznacz[it]+ "--")

    ps.ticks_changes(ax)

#reading model output from text file and converting to an array
def reading_modeloutput(filename):
    f = h5py.File(filename, "r")
    h = np.array(f["h"])
    qx = np.array(f["qx"])
    qy = np.array(f["qy"])
    return h, qx, qy 

#plotting together analytic solution and model output 
def analytic_model_fig(ax, x_range, y_range, h_m, v_m, time=1):
    lamb = eq.d2_lambda_evol(time)
    h_a = eq.d2_height(lamb, x_range, y_range)
    v_a = eq.d2_velocity(lamb, x_range, y_range)

    ax.plot(x_range, eq.d2_initial(x_range, y_range), 'k', x_range, h_a, 'b',
            x_range, h_m, "r")
    ax.plot(x_range, 0*x_range, "k-", x_range, v_a, 'b--',
            x_range, v_m, "r--")

    #ax.set_ylim(-2,2)
    ps.ticks_changes(ax)


# time_l - list of time levels for analytic solutions
# time - model time level used for comparison with analytic solution
# dt -  model time step
# x_shift - shift between initial cond. in model and for analytic solution
def main(dir, casename_l, x_shift=8, time_l=[0,3], time=3, dt=0.01):
    plt.figure(1, figsize = (6,3))
    ax = plt.subplot(1,1,1)
    #plotting analytic solution
    analytic_fig(ax, time_l)
    plt.savefig("papier_shallowwater_2d_analytic.pdf")
    #plt.show()
    #plotting comparison between analytic solution and model results for various options
    for ic, casename in enumerate(casename_l):
        plt.figure(ic+1, figsize = (6,3))
        ax = plt.subplot(1,1,1)

        print "plotting for " + casename + ", t = " + str(time)
        #model variables TODO: time
        h_m, px_m, py_m = reading_modeloutput(dir+"spreading_drop_2d_" + casename +
                                              ".out/timestep0000000" + str(int(time/dt))
                                              + '.h5')
        # calculate velocity (momentum/height) only in the droplet region.
        #calculating velocity from momentum, only for the droplet area 
        v_m = np.where(h_m > 0,  py_m/h_m, 0)
        print "where with h_m = 0 !!"
                
        # choosing a plane of a cross section TODO: should be 160?
        # TODO: x_range/y_range should be calculated from hdf file!!
        print "TODO: x_range/y_range should be calculated from hdf file!!"
        analytic_model_fig(ax,
                           np.linspace(-8,8,h_m.shape[0]), np.zeros(h_m[0].shape[0]),
                           h_m[159], v_m[159], time)
        
        #ax.annotate(str(casename), xy=(0.01, 0.97), xycoords='axes fraction',
        #            fontsize=12, horizontalalignment='left', verticalalignment='top')
        
        plt.savefig("papier_shallowwater_2d_"+str(casename)+".pdf")
        #plt.show()

main("./", sys.argv[1:])

    
    
    
