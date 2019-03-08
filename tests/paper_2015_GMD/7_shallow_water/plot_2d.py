import numpy as np
import h5py
import sys

import matplotlib
matplotlib.use('Pdf')
import matplotlib.pyplot as plt

import plot_settings as ps
import analytic_eq as eq

# plotting analytic solutions for height and velocity
def analytic_fig(ax, x_lim, time_l=[0,1,2,3], nxy=320):
    x_range = np.linspace(-x_lim, x_lim, nxy)
    y_range = np.zeros(nxy)
    oznacz = ['k', 'b', 'c', 'y', 'g', 'm', 'r']
    for it, time in enumerate(time_l):
        lamb = eq.d2_lambda_evol(time)
        h = eq.d2_height(lamb, x_range, y_range)
        v = eq.d2_velocity(lamb, x_range, y_range)
        ax.plot(x_range, h, oznacz[it])
        ax.plot(x_range, v, oznacz[it]+ "--")
    ps.ticks_changes(ax)

# reading model output from text file and converting to an array
def reading_modeloutput(dir, time):
    dir_model = {}
    f_crd = h5py.File(dir+ "/const.h5", "r")
    time_model = np.array(f_crd["T"])
    assert(time in time_model),"time level not in model output"
    dt = round(f_crd["advection"].attrs["dt"][0], 4)
    dir_model["dx"] = f_crd["advection"].attrs["di"][0]
    f_out = h5py.File(dir+"/timestep0000000" + str(int(time/dt))+ '.h5', "r")
    dir_model["h"] = np.array(f_out["h"])
    dir_model["qx"] = np.array(f_out["qx"])
    dir_model["qy"] = np.array(f_out["qy"])
    return dir_model

# plotting together analytic solution and model output
def analytic_model_fig(ax, var_md, time=1):
    x_range = var_md["x_range"]
    y_range = 0 * var_md["x_range"]
    ind_cs = int(var_md["h"].shape[1]/2)
    lamb = eq.d2_lambda_evol(time)
    h_a = eq.d2_height(lamb, x_range, y_range)
    v_a = eq.d2_velocity(lamb, x_range, y_range)
    ax.plot(x_range, eq.d2_initial(x_range, y_range), 'k', x_range, h_a, 'b',
            x_range, var_md["h"][:,ind_cs], "r")
    ax.plot(x_range, 0*x_range, "k-", x_range, v_a, 'b--',
            x_range, var_md["vx"][:,ind_cs], "r--")
    ps.ticks_changes(ax)

# time_an - list of time levels for analytical solutions
# time - model time level used for comparison with analytic solution
# xy_lim - limit used for calculating x_range for analytical solution
def main(dir, casename_l, xy_lim=8, time_l=[0,3], time=3):
    plt.figure(1, figsize = (6,3))
    ax = plt.subplot(1,1,1)
    # plotting analytic solution
    analytic_fig(ax, xy_lim, time_l)
    plt.savefig("2d_analytic.pdf")
    # plotting comparison between analytic solution and model results for various options
    for ic, casename in enumerate(casename_l):
        plt.figure(ic+2, figsize = (6,3))
        ax = plt.subplot(1,1,1)
        print("plotting for " + casename + ", t = " + str(time))
        var_model = reading_modeloutput(dir + casename, time)
        # calculate velocity (momentum/height) only in the droplet region.
        var_model["vx"] = np.where(var_model["h"]>0, var_model["qx"]/var_model["h"], 0)
        # calculating model output coord.
        var_model["x_range"] = np.arange(-xy_lim+var_model["dx"]/2., xy_lim, var_model["dx"])
        assert(var_model["x_range"].shape[0] == var_model["h"].shape[0]), "domain size differs from model output shape"

        analytic_model_fig(ax, var_model, time)
        plt.savefig(str(casename)+".pdf")

main("./", sys.argv[1:])
