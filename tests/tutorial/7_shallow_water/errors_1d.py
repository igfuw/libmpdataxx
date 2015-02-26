import sys
sys.path.append("../analytical/")
import numpy as np
import h5py
import analytic_eq as eq

#reading model output and saving as numpy arrays   
def reading_modeloutput(filename):
    f  = h5py.File(filename, "r")
    h  = np.array(f["h"]).transpose()
    qx = np.array(f["qx"]).transpose()
    return h, qx

def errors(dir, dt, dx, time_l, x_lim):
    x_range = np.arange(-x_lim, x_lim, dx)

    for it, time in enumerate(time_l):

        # reading the model output
        time_str = '%0*d' % (10, int(time/dt))
        print dir + "/timestep" + time_str + '.h5'
        h_m, px_m = reading_modeloutput(dir + "/timestep" + time_str + '.h5')

        # calculating the analytical solution 
        lamb = eq.d1_lambda_evol(time)
        h_an = eq.d1_height(lamb, x_range)

        # calculating the errors of the drop depth, eq. 25 and 26 from the paper
        h_diff    = h_m - h_an
        points_nr = h_m.shape[0]
        delh_inf  = abs(h_diff).max() 
        delh_2    = 1./time * ((h_diff**2).sum() / points_nr )**0.5 

        # outputing general info
        if time == 1:
            file = open(dir + "_stats.txt", "w")
            file.write( "dx                                 = " + str(dx)                       + "\n")
            file.write( "dt                                 = " + str(dt)                       + "\n")
            file.write( "number of points in the domain     = " + str(points_nr)                + "\n")
            file.write( "L_inf                              = max|h_m-h_an|"                    + "\n")
            file.write( "L_2                                = sqrt(sum(h_m-h_an)^2 / N) / time" + "\n" + "\n")

        # outputting error statistics
        file.write( "time                               = " + str(time)                         + "\n")
        file.write( "max(h_an)                          = " + str(round(h_an.max(), 8))         + "\n")
        file.write( "max(h_num)                         = " + str(round(h_m.max(), 8))          + "\n")
        file.write( "L_inf                              = " + str(round(delh_inf, 8))           + "\n")
        file.write( "L_2                                = " + str(round(delh_2, 8))             + "\n" + "\n")
        if time == 3:
            file.close()
        
# printing errors at different time steps
def evolution_test(dir, dt, dx, time_l=[1,2,3], x_lim=8):
    errors(dir, dt, dx, time_l, x_lim)
    
def main(dir, casename_l):
    for casename in casename_l:
        print casename
        evolution_test(dir + str(casename), dt=0.01, dx=0.05) #TODO: read it from the h5 file
        
main("./", sys.argv[1:])
