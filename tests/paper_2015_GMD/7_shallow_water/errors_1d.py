import sys
sys.path.append("../analytical/")
import numpy as np
import h5py
import analytic_eq as eq

#reading model output and saving as numpy arrays
def reading_modeloutput(dir, time):
    dir_model = {}
    f_crd = h5py.File(dir+ "/const.h5", "r")
    time_model = np.array(f_crd["T"])
    assert(time in time_model),"time level not in model output"
    dir_model["dt"] = round(f_crd["advection"].attrs["dt"][0], 4)
    dir_model["dx"] = f_crd["advection"].attrs["di"][0]
    f_out = h5py.File(dir+"/timestep0000000" + str(int(time/dir_model["dt"]))+ '.h5', "r")
    dir_model["h"] = np.array(f_out["h"])
    dir_model["qx"] = np.array(f_out["qx"])
    return dir_model

def errors(dir, time_l, x_lim):
    for it, time in enumerate(time_l):
        # reading the model output
        var_model = reading_modeloutput(dir, time)

        var_model["x_range"] = np.arange(-x_lim+var_model["dx"]/2., x_lim, var_model["dx"])
        assert(var_model["x_range"].shape == var_model["h"].shape), "domain size differs from model output shape"

        # calculating the analytical solution
        lamb = eq.d1_lambda_evol(time)
        h_an = eq.d1_height(lamb, var_model["x_range"])

        # calculating the errors of the drop depth, eq. 25 and 26 from the paper
        h_diff    = var_model["h"] - h_an
        points_nr = var_model["h"].shape[0]
        delh_inf  = abs(h_diff).max()
        delh_2    = 1./time * ((h_diff**2).sum() / points_nr )**0.5

        # outputing general info
        if time == time_l[0]:
            file = open(dir + "_stats.txt", "w")
            fstring = (
                       "dx                                 = {dx:.4f}\n"
                       "dt                                 = {dt:.4f}\n"
                       "number of points in the domain     = {npoints}\n"
                       "L_inf                              = max|h_m-h_an|\n"
                       "L_2                                = sqrt(sum(h_m-h_an)^2 / N) / time\n\n"
                      )
            file.write(fstring.format(dx = var_model["dx"], dt = var_model["dt"], npoints = points_nr))

        # outputting error statistics
        fstring = (
                   "time                               = {time:.4f}\n"
                   "max(h_an)                          = {max_h_an:.8f}\n"
                   "max(h_num)                         = {max_h_num:.8f}\n"
                   "L_inf                              = {L_inf:.8f}\n"
                   "L_2                                = {L_2:.8f}\n"
                   "max(px_num)                        = {max_px_num:.8f}\n\n"
                  )
        file.write(fstring.format(
                                  time = time,
                                  max_h_an = h_an.max(),
                                  max_h_num = var_model["h"].max(),
                                  L_inf = delh_inf,
                                  L_2 = delh_2,
                                  max_px_num = var_model["qx"].max()
                                 )
                  )
    file.close()

# printing errors at different time steps
def evolution_test(dir, time_l=[1,2,3], x_lim=8):
    errors(dir, time_l, x_lim)

def main(dir, casename_l):
    for casename in casename_l:
        print(casename)
        evolution_test(dir + str(casename)) #TODO: read it from the h5 file

main("./", sys.argv[1:])
