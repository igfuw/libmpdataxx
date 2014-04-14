from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import h5py
import plot_settings as ps

def reading_modeloutput(filename):
    f = h5py.File(filename, "r")
    h = np.array(f["h"])
    qx = np.array(f["qx"])
    qy = np.array(f["qy"])
    return h, qx, qy 



def plotting_2D(X, Y, Z):
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.plot_surface(X, Y, Z, rstride=10, cstride=10, alpha=0.2)
    cset = ax.contourf(X, Y, Z, zdir='z', offset=-0.1, cmap=cm.Blues)
    fig.colorbar(cset) #, fraction=0.05)
    ps.ticks_changes(ax)


    ax.set_xlabel('x')
    #ax.set_xlim(-40, 40)
    ax.set_ylabel('y')
    #ax.set_ylim(-40, 40)
    ax.set_zlabel('h')
    ax.set_zlim(-0.1, 0.1)

    plt.savefig("plot2D_it=3.pdf")
    plt.show()


#TODO napisac uwage o dx, dy
def main(filename):
    x_range = np.linspace(-8,8,320)
    y_range = np.linspace(-8,8,320)
    X, Y = np.meshgrid(x_range, y_range)
    h, qx, qy = reading_modeloutput(filename)
    plotting_2D(X, Y, h)

main("./spreading_drop_2d_fct+iga.out/timestep0000000300.h5")
