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
    f_crd = h5py.File(dir+ "/coord.h5", "r")
    time_model = np.array(f_crd["T"])
    assert(time in time_model),"time level not in model output"
    dt = round(f_crd["T"].attrs["dt"], 4)
    dir_model["X"] = np.array(f_crd["X"])
    dir_model["Y"] = np.array(f_crd["Y"])
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


def main(dir, time):
    var_model = reading_modeloutput(dir, time)
    # model coord.
    x_model_shift = 0.5 * var_model["X"][:,0].max()
    var_model["x_range"] = 0.5*(var_model["X"][1:,0]+var_model["X"][:-1,0]) - x_model_shift
    y_model_shift = 0.5 * var_model["Y"][0,:].max()
    var_model["y_range"] = 0.5*(var_model["Y"][0,1:]+var_model["Y"][0,:-1]) - y_model_shift
    X, Y = np.meshgrid(var_model["x_range"], var_model["y_range"])

    plotting_3D(X, Y, var_model["h"])

main(dir="./2d_fct_iga", time=3)
