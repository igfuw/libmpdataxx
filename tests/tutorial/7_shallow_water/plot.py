import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib.ticker as ticker
import numpy as np
import sys
import plot_settings as ps
import analytic_eq as eq

# plotting analytic solutions for height and velocity 
def analytic_fig(ax, time_l = [0,1,2,3], nx=320):
    x_range = np.linspace(-8,8, nx)
    oznacz = ['k', 'b', 'c', 'y', 'g', 'm', 'r']
    y0 = eq.d1_initial(x_range)
    for it, time in enumerate(time_l):
        lamb = eq.d1_lambda_evol(time)
        h = eq.d1_height(lamb, x_range)
        v = eq.d1_velocity(lamb, x_range)
        ax.plot(x_range, h, oznacz[it])
        ax.plot(x_range, v, oznacz[it]+ "--")
    ps.ticks_changes(ax)

#reading model output from text file and converting to an array
def reading_modeloutput(filename):
    f = open(filename)
    list_nt = []
    for line in f:
        lst = line.split('\t')
        list_nt.append(lst[:-1])
    return np.array(list_nt, dtype=float) 

#plotting together analytic solution and model output 
def analytic_model_fig(ax, x_range, h_m, v_m, t_m, it):
    lamb = eq.d1_lambda_evol(t_m[it,0])
    h_a = eq.d1_height(lamb, x_range)
    v_a = eq.d1_velocity(lamb, x_range)

    ax.plot(x_range, eq.d1_initial(x_range), 'k', x_range, h_a, 'b',
            x_range, h_m[it], "r")
    ax.plot(x_range, 0*x_range, "k-", x_range, v_a, 'b--',
            x_range, v_m[it], "r--")

    ax.set_ylim(-2,2)
    ps.ticks_changes(ax)

# time_l - list of time levels for analytic solutions
# it - model time level used for comparison with analytic solution
# x_shift - shift between initial cond. in model and for analytic solution
def main(dir, casename_l, x_shift=8, time_l=[0,3], it=300):
    plt.figure(1, figsize = (6,3))
    ax = plt.subplot(1,1,1)
    #plotting analytic solution
    print "plotting analytic solution"
    analytic_fig(ax, time_l)
    plt.savefig("papier_shallowwater_1d_analytic.pdf")
    #plt.show()
    #plotting comparison between analytic solution and model results for various options
    for ic, casename in enumerate(casename_l): 
        plt.figure(ic+1, figsize = (6,3))
        ax = plt.subplot(1,1,1)
        h_m = reading_modeloutput(dir+casename+".h")
        p_m = reading_modeloutput(dir+casename+".q")
        #calculating velocity from momentum, only for the droplet area 
        v_m = np.where(h_m > 0,  p_m/h_m, 0)
        print "where with h_m = 0 !!"
        x_m = reading_modeloutput(dir+casename+".x") - x_shift
        t_m = reading_modeloutput(dir+casename+".t")
        x_range  = x_m[0]

        print "plotting " + str(casename) + ", t = " + str(t_m[it,0])        
        analytic_model_fig(ax, x_range, h_m, v_m, t_m, it)
        # annotating figures
        # ax.annotate(str(casename), xy=(0.01, 0.97), xycoords='axes fraction',
        #             fontsize=12, horizontalalignment='left', verticalalignment='top')
        

        plt.savefig("papier_shallowwater_1d_"+str(casename)+".pdf")
        #plt.show()


main("./", sys.argv[1:])

    
    
    
