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
    dir_model["dy"] = f_crd["advection"].attrs["dj"][0]
    f_out = h5py.File(dir+"/timestep0000000" + str(int(time/dir_model["dt"]))+ '.h5', "r")
    dir_model["h"] = np.array(f_out["h"])
    dir_model["qx"] = np.array(f_out["qx"])
    dir_model["qy"] = np.array(f_out["qy"])
    return dir_model

def errors(dir, time_l, xy_lim):
    for time in time_l:
        # reading the model output
        var_model = reading_modeloutput(dir, time)

        var_model["x_range"] = np.arange(-xy_lim+var_model["dx"]/2., xy_lim, var_model["dx"])
        assert(var_model["x_range"].shape[0] == var_model["h"].shape[0]), "domain size differs from model output shape"
        var_model["y_range"] = np.arange(-xy_lim+var_model["dy"]/2., xy_lim, var_model["dx"])
        assert(var_model["y_range"].shape[0] == var_model["h"].shape[1]), "domainsize differs from model output shape"

        # calculating the analytical solution
        lamb = eq.d2_lambda_evol(time)
        h_an = eq.d2_height_plane(lamb, var_model["x_range"], var_model["y_range"])

        # calculating the errors of the drop depth, eq. 25 and 26 from the paper
        h_diff    = var_model["h"] - h_an
        points_nr = var_model["h"].shape[0] * var_model["h"].shape[1]
        delh_inf  = abs(h_diff).max()
        delh_2    = 1./time * ((h_diff**2).sum() / points_nr )**0.5

        # outputing general info
        if time == time_l[0]:
            file = open(dir + "_stats.txt", "w")
            fstring = (
                       "dx                                 = {dx:.4f}\n"
                       "dy                                 = {dy:.4f}\n"
                       "dt                                 = {dt:.4f}\n"
                       "number of points in the domain     = {npoints}\n"
                       "L_inf                              = max|h_m-h_an|\n"
                       "L_2                                = sqrt(sum(h_m-h_an)^2 / N) / time\n\n"
                      )
            file.write(fstring.format(
                                      dx = var_model["dx"],
                                      dy = var_model["dy"],
                                      dt = var_model["dt"],
                                      npoints = points_nr
                                     )
                      )

        # outputting error statistics
        fstring = (
                   "time                               = {time:.4f}\n"
                   "max(h_an)                          = {max_h_an:.7f}\n"
                   "max(h_num)                         = {max_h_num:.7f}\n"
                   "L_inf                              = {L_inf:.7f}\n"
                   "L_2                                = {L_2:.7f}\n"
                   "max(px_num)                        = {max_px_num:.7f}\n"
                   "max(py_num)                        = {max_py_num:.7f}\n\n"
                  )
        file.write(fstring.format(
                                  time = time,
                                  max_h_an = h_an.max(),
                                  max_h_num = var_model["h"].max(),
                                  L_inf = delh_inf,
                                  L_2 = delh_2,
                                  max_px_num = var_model["qx"].max(),
                                  max_py_num = var_model["qy"].max()
                                 )
                  )
    file.close()

# printing errors at different time steps
def evolution_test(dir, time_l=[1,2,3], xy_lim=8):
    errors(dir, time_l, xy_lim)

def main(dir, casename_l):
    for casename in casename_l:
        print(casename)
        evolution_test(dir + str(casename))

main("./", sys.argv[1:])
