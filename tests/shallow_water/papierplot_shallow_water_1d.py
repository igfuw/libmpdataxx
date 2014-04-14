from scipy.optimize import fsolve
import math
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib.ticker as ticker
import numpy as np
import sys


#eq. 5.6 (Schar&Smolarkiewicz, 1996) 
def lambda_eq(x, *time):
    return 0.5 * ((x*(x-1.))**0.5 + math.log((x-1.)**0.5 + x**0.5)) - time

# finging roots of lambda_eq 
def lambda_evol(time, x0=1):
    return fsolve(lambda_eq, x0, args=time)

#eq. 5.4 (Schar&Smolarkiewicz, 1996)
def height(lamb, x):
    return np.where(x**2 <= lamb**2, lamb**-1 * (1. - (x/lamb)**2), 0)

#eq. 5.5 (Schar&Smolarkiewicz, 1996)
def velocity(lamb, x):
    lamb_t = 2 * (1 - lamb**-1)**0.5
    return np.where(x**2 <= lamb**2, x * lamb_t / lamb, 0)

#eq. 5.3 (Schar&Smolarkiewicz, 1996)
def initial(x):
    return np.where(x**2<=1, 1-x**2, 0)

# plotting analytic solutions for height and velocity 
def analytic_fig(ax,time_l = [0,1,2,3], nx=320):
    x_range = np.linspace(-8,8, nx)
    oznacz = ['k', 'b', 'c', 'y', 'g', 'm', 'r']
    y0 = initial(x_range)

    for it, time in enumerate(time_l):
        lamb = lambda_evol(time)
        h = height(lamb, x_range)
        v = velocity(lamb, x_range)
        ax.plot(x_range, h, oznacz[it])
        ax.plot(x_range, v, oznacz[it]+ "--")

    # removing some ticks' labels
    for i, tick in enumerate(ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks()):
        if i % 2 != 0:
            tick.label1On = False                                

    # changing ticks' size  
    for item in plt.xticks()[1] + plt.yticks()[1]:
        item.set_fontsize(10)



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
    lamb = lambda_evol(t_m[it,0])
    h_a = height(lamb, x_range)
    v_a = velocity(lamb, x_range)

    ax.plot(x_range, initial(x_range), 'k', x_range, h_a, 'b',
            x_range, h_m[it], "r")
    ax.plot(x_range, 0*x_range, "k--", x_range, v_a, 'b--',
            x_range, v_m[it], "r--")

    ax.set_ylim(-2,2)

    # removing some ticks' labels
    for i, tick in enumerate(ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks()):
        if i % 2 != 0:
            tick.label1On = False                                

    # changing ticks' size  
    for item in plt.xticks()[1] + plt.yticks()[1]:
        item.set_fontsize(10)

# time_l - list of time levels for analytic solutions
# it - model time level used for comparison with analytic solution
# x_shift - shift between initial cond. in model and for analytic solution
def main(dir, casename_l, x_shift=8, time_l=[0,3], it=300):
    plt.figure(1, figsize = (6,8))
    ax = plt.subplot(len(casename_l)+1,1,1)
    #plotting analytic solution
    print "plotting analytic solution"
    analytic_fig(ax, time_l)
    #plotting comparison between analytic solution and model results for various options
    for ic, casename in enumerate(casename_l): 
        h_m = reading_modeloutput(dir+casename+".h")
        p_m = reading_modeloutput(dir+casename+".q")
        #calculating velocity from momentum, only for the droplet area 
        v_m = np.where(h_m > 0,  p_m/h_m, 0)
        print "where with h_m = 0 !!"
        x_m = reading_modeloutput(dir+casename+".x") - x_shift
        t_m = reading_modeloutput(dir+casename+".t")
        x_range  = x_m[0]

        print "plotting " + str(casename) + ", t = " + str(t_m[it,0])        
        ax = plt.subplot(len(casename_l)+1,1,ic+2)
        analytic_model_fig(ax, x_range, h_m, v_m, t_m, it)
        # annotating figures
        #ax.annotate(str(casename), xy=(0.01, 0.97), xycoords='axes fraction',
        #            fontsize=12, horizontalalignment='left', verticalalignment='top')
        
    plt.savefig("papier_shallowwater_1d.pdf")
    plt.show()

main("./", sys.argv[1:])

    
    
    
