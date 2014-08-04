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
def analytic_fig(ax, time_l = [0.,3., 7.], x_range = np.linspace(-8,8,320),
                              y_range = np.zeros(320)):
    oznacz = ['k', 'b', 'c', 'y', 'g', 'm', 'r']
    y0 = eq.d2_initial(x_range, y_range)

    for it, time in enumerate(time_l):
        print "time w fig", time
        if time == 0.:
            lamb_ar = [3., 0., 1., 0.]
        else:
            lamb_ar = eq.d2_el_lamb_lamb_t_evol(time, lamb_x0=3., lamb_y0=1.)#TODO
        h_x = eq.d2_el_height(lamb_ar[0], lamb_ar[2], x_range, y_range)
        #TODO zrobic to lepiej, zeby sie nie mylilo x z y
        h_y = eq.d2_el_height(lamb_ar[0], lamb_ar[2], y_range, x_range)
        v_x = eq.d2_el_velocity(lamb_ar[0], lamb_ar[1], lamb_ar[2], x_range, y_range)
        v_y = eq.d2_el_velocity(lamb_ar[2], lamb_ar[3], lamb_ar[0], x_range, y_range)
        import pdb
        #pdb.set_trace()
        ax.plot(x_range, h_x, oznacz[it])
        ax.plot(x_range, v_x, oznacz[it]+ "--")
        ax.set_ylim(-1.,1.)

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
    print "time analytic_model", time
    lamb_ar = eq.d2_el_lamb_lamb_t_evol(time, lamb_x0=3., lamb_y0=1.)#TODO
    h_x = eq.d2_el_height(lamb_ar[2], lamb_ar[0], x_range, y_range) #TODO: zmienia sie
    h_x0 = eq.d2_el_height(1., 3., x_range, y_range) # TODO: zmienia sie
    #TODO zrobic to lepiej, zeby sie nie mylilo x z y
    h_y = eq.d2_el_height(lamb_ar[0], lamb_ar[2], y_range, x_range)
    h_y0 = eq.d2_el_height(3., 1., y_range, x_range)
    v_x = eq.d2_el_velocity(lamb_ar[0], lamb_ar[1], lamb_ar[2], x_range, y_range)
    v_y = eq.d2_el_velocity(lamb_ar[2], lamb_ar[3], lamb_ar[0], x_range, y_range)

    ax.plot(x_range, h_x0, 'k', x_range, h_x, 'b',
            x_range, h_m, "r")
    ax.plot(x_range, 0*x_range, "k--", x_range, v_y, 'b--', #TODO: zmienia sie
            x_range, v_m, "r--")

    ax.set_ylim(-1., 1.)
    ps.ticks_changes(ax)


# time_l - list of time levels for analytic solutions
# time - model time level used for comparison with analytic solution
# dt -  model time step
# x_shift - shift between initial cond. in model and for analytic solution
def main(dir, casename_l, x_shift=8, time_l=[0,3], time=1, dt=0.01):
    plt.figure(1, figsize = (6,8))
    ax = plt.subplot(len(casename_l)+1,1,1)
    #plotting analytic solution
    analytic_fig(ax)
    #plotting comparison between analytic solution and model results for various options
    for ic, casename in enumerate(casename_l):
        print "plotting for " + casename + ", t = " + str(time)
        #model variables TODO: time
        print "filename", (dir+"spreading_drop_2delipsa_" + casename +
                                              ".out/timestep0000000" + str(int(time/dt))
                                              + '.h5')
        
        h_m, px_m, py_m = reading_modeloutput(dir+"spreading_drop_2delipsa_" + casename +
                                              ".out/timestep0000000" + str(int(time/dt))
                                              + '.h5')
        # calculate velocity (momentum/height) only in the droplet region.
        #calculating velocity from momentum, only for the droplet area 
        v_m = np.where(h_m > 0.,  py_m/h_m, 0)
        print "where with h_m = 0 !!"
        import pdb
        pdb.set_trace()
        ax = plt.subplot(len(casename_l)+1,1,ic+2)
        
        # choosing a plane of a cross section TODO: should be 160?
        # TODO: x_range/y_range should be calculated from hdf file!!
        print "TODO: x_range/y_range should be calculated from hdf file!!"
        analytic_model_fig(ax,
                           np.linspace(-8,8,h_m.shape[0]), np.zeros(h_m[0].shape[0]),
                           h_m[159,:], v_m[159,:], time)
        
        ax.annotate(str(casename), xy=(0.01, 0.97), xycoords='axes fraction',
                    fontsize=12, horizontalalignment='left', verticalalignment='top')
        
    plt.savefig("papier_shallowwater_2delipsa_x=0_it=1.pdf")
    plt.show()

main("./", sys.argv[1:])

    
    
    
