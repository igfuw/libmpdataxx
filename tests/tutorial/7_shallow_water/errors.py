import sys
sys.path.append("../analytical/")
import numpy as np
import h5py
import analytic_eq as eq

#reading model output and saving as numpy arrays   
def reading_modeloutput(filename):
    f = h5py.File(filename, "r")
    h = np.array(f["h"]).transpose()
    qx = np.array(f["qx"]).transpose()
    qy = np.array(f["qy"]).transpose()
    return h, qx, qy


def errors(dir, dt, dx, time_l, xy_lim):
    x_range = y_range = np.arange(-xy_lim, xy_lim, dx)
    for it, time in enumerate(time_l):
        print "\n", "TIME t = " + str(time), "dx, dt", dx, dt
        time_str = '%0*d' % (10, int(time/dt))
        # reading the model output
        h_m, px_m, py_m = reading_modeloutput(dir +  "/timestep" + time_str + '.h5')
        # calculating the analytical solution 
        lamb = eq.d2_lambda_evol(time)
        h_an = eq.d2_height_plane(lamb, x_range, y_range)
                        
        # calculating the errors of the drop depth, eq. 25 and 26 from the paper
        h_diff = h_m - h_an
        points_nr = h_m.shape[0]*h_m.shape[1]
        print "number of points in the domain, max(h_an), max(h_num)", points_nr, h_an.max(), h_m.max()
        delh_inf = abs(h_diff).max() 
        delh_2 = 1./time * ((h_diff**2).sum() / points_nr )**0.5 

        print "L_inf = max|h_m-h_an| = ", delh_inf 
        print "L_2 = sqrt(sum(h_m-h_an)^2 / N) / time = ", delh_2, "\n"
        
# comparing errors for reference simulation at different time steps
def evolution_test(dir, dt, dx, time_l=[1,2,3], xy_lim=8):
    errors(dir, dt, dx, time_l, xy_lim)
    

def main(dir, casename_l):
    for casename in casename_l:
        print "****** errrors for " + casename
        evolution_test(dir + "spreading_drop_2d_" + str(casename) + ".out",
                       dt=0.01, dx=0.05)
        
main("./", sys.argv[1:])
