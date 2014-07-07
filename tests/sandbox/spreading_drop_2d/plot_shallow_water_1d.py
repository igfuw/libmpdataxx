from scipy.optimize import fsolve
import math
from pylab import *
from matplotlib.font_manager import FontProperties
import matplotlib.ticker as ticker
import numpy as np


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
def height_veloc_analytic_fig(time_l = [0,1,2,3]):
    x_range = np.linspace(-8,8,100)
    oznacz = ['c', 'm', 'k', 'y', 'b', 'g', 'r']
    y0 = initial(x_range)

    figure(1, figsize = (8,4))
    ax = subplot(1,1,1)

    for it in time_l:
        lamb = lambda_evol(it)
        h = height(lamb, x_range)
        v = velocity(lamb, x_range)
        ax.plot(x_range, h, oznacz[it])
        ax.plot(x_range, v, oznacz[it]+ "--")

    savefig("height_velocity_analytic.pdf")
    show()

#reading model output from text file and converting to an array
def reading_modeloutput(filename):
    f = open(filename)
    list_nt = []
    for line in f:
        lst = line.split('\t')
        list_nt.append(lst[:-1])
    return np.array(list_nt, dtype=float) 

#plotting together analytic solution and model output 
def analytic_model_fig(x_range, h_m, v_m, t_m, casename):
    for it in range(t_m.shape[0]): #TODO moze ladniej?
        lamb = lambda_evol(t_m[it,0])
        h_a = height(lamb, x_range)
        v_a = velocity(lamb, x_range)

        figure(1, figsize = (4,3))
        ax = subplot(1,1,1)

        ax.plot(x_range, initial(x_range), 'k', x_range, h_a, 'b',
                x_range, h_m[it], "r")
        ax.plot(x_range, 0*x_range, "k--", x_range, v_a, 'b--',
                x_range, v_m[it], "r--")

        ax.set_title("t = "+str(t_m[it,0]), fontsize=12)

        ax.set_ylim(-2,2)

        # removing some ticks labels
        for i, tick in enumerate(ax.xaxis.get_major_ticks()+ax.yaxis.get_major_ticks()):
            if i % 2 != 0:
                tick.label1On = False                                
        
        # changing size of all ticks  
        for item in xticks()[1] + yticks()[1]:
            item.set_fontsize(10)

        savefig("compar"+casename+"_height_velocity_t="+str(t_m[it,0])+".pdf")
        show()
    
    

def main(dir, casename, x_shift=8, time_l=[1,2,3]):
    #plotting anlytic solutions for height and velocity
    height_veloc_analytic_fig()

    #model variables
    h_m_nt = reading_modeloutput(dir+"spreading_drop_1d"+casename+".out/out.h")
    p_m_nt = reading_modeloutput(dir+"spreading_drop_1d"+casename+".out/out.q")
    v_m_nt = np.where(p_m_nt != 0.,  p_m_nt/h_m_nt, 0)
    print "where with p != 0 !!"
    x_m = reading_modeloutput(dir+"spreading_drop_1d"+casename+".out/out.x") - x_shift
    t_m = reading_modeloutput(dir+"spreading_drop_1d"+casename+".out/out.t")

    #finding indices for time levels in time_l
    #TODO: informacja, ze nie ma danej chwili czasu
    ind_time = np.in1d(t_m,time_l)

   
    #plotting analytic solutions and model outputs for time from time_l
    analytic_model_fig(x_m[0], h_m_nt[ind_time], v_m_nt[ind_time], t_m[ind_time],
                       casename)

main("../../build/tests/shallow_water/", "_fct_iga")
    
    
    
