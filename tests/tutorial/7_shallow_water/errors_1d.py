import sys
sys.path.append("../analytical/")
import numpy as np
import h5py
import analytic_eq as eq

#reading model output from text file and converting to an array                
def reading_modeloutput(filename):
    f = open(filename)
    list_nt = []
    for line in f:
        lst = line.split('\t')
        list_nt.append(lst[:-1])
    return np.array(list_nt, dtype=float)



def errors(dir, casename, dt_output, dx, time_l, x_lim):
    x_range = np.arange(-x_lim, x_lim, dx)
    # reading the model output (all time levels) 
    h_m_nt = reading_modeloutput(dir+casename+".h")
    p_m_nt = reading_modeloutput(dir+casename+".q")
    t_m_nt = reading_modeloutput(dir+casename+".t")

    for it, time in enumerate(time_l):
        it_output = int(time/dt_output)
        print "\n", "TIME t = " + str(time)
        assert (time == t_m_nt[it_output][0]), "something wrog with time levels"

        # calculating the analytical solution 
        lamb = eq.d1_lambda_evol(time)
        h_an = eq.d1_height(lamb, x_range)
        # calculating the errors of the drop depth, eq. 25 and 26 from the paper
        h_diff = h_m_nt[it_output] - h_an
        points_nr = h_m_nt[it_output].shape[0]
        print "number of points in the domain, max(h_an), max(h_num)", points_nr, h_an.max(), h_m_nt[it_output].max()
        delh_inf = abs(h_diff).max() 
        delh_2 = 1./time * ((h_diff**2).sum() / points_nr )**0.5 

        print "L_inf = max|h_m-h_an| = ", delh_inf 
        print "L_2 = sqrt(sum(h_m-h_an)^2 / N) / time = ", delh_2, "\n"
        
# comparing errors for reference simulation at different time steps
def evolution_test(dir, casename, dt_output, dx, time_l=[1,2,3], x_lim=8):
    errors(dir, casename, dt_output, dx, time_l, x_lim)
    

def main(dir, casename_l):
    for casename in casename_l:
        print "****** errrors for " + casename
        evolution_test(dir, casename,
                       dt_output=0.01, dx=0.05)
        
main("./", sys.argv[1:])
