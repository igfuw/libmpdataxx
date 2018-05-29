import matplotlib
matplotlib.use('Pdf')
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import h5py
import plot_settings as ps


# reading model output from text file and converting to an array
def reading_modeloutput(dir, time):
    dir_model = {}
    f_crd = h5py.File(dir+ "/const.h5", "r")
    time_model = np.array(f_crd["T"])
    assert(time in time_model),"time level not in model output"
    dt = round(f_crd["advection"].attrs["dt"][0], 4)
    dir_model["dx"] = f_crd["advection"].attrs["di"][0]
    dir_model["dy"] = f_crd["advection"].attrs["dj"][0]
    f_out = h5py.File(dir+"/timestep0000000" + str(int(time/dt))+ '.h5', "r")
    dir_model["h"] = np.array(f_out["h"])
    return dir_model

def plotting_3D(X, Y, Z):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, Z, rstride=10, cstride=10, alpha=0.2)
    cset = ax.contourf(X, Y, Z, zdir='z', offset=-0.1, cmap=cm.Blues)
    fig.colorbar(cset)
    ps.ticks_changes(ax)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('h')
    ax.set_zlim(-0.1, 0.1)
    plt.savefig("2d_fct_iga_view.pdf")

def main(dir, time, xy_lim=8):
    var_model = reading_modeloutput(dir, time)
    # calculating model output coord.
    var_model["x_range"] = np.arange(-xy_lim+var_model["dx"]/2., xy_lim, var_model["dx"])
    var_model["y_range"] = np.arange(-xy_lim+var_model["dy"]/2., xy_lim, var_model["dy"])
    assert(var_model["x_range"].shape[0] == var_model["h"].shape[0]), "domain size differs from model output shape"
    assert(var_model["y_range"].shape[0] == var_model["h"].shape[1]), "domain size differs from model output shape"

    X, Y = np.meshgrid(var_model["x_range"], var_model["y_range"])
    plotting_3D(X, Y, var_model["h"])

main(dir="./2d_fct_iga", time=3)
