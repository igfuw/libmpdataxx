## @file
# @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
# @copyright University of Warsaw
# @date Januar 2012
# @section LICENSE
# GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)

# 1d isolines test from pks & wwg 1990

import numpy as np # arrays
import subprocess # shell calls
import os # unlink()
import sys # exit, args
import math # sqrt
import matplotlib.pyplot as plt # plots
from matplotlib.mlab import griddata # griddata

dx, cour, err = np.loadtxt(sys.argv[1], unpack=True)

theta=np.zeros(dx.shape[0])
r=np.zeros(dx.shape[0])
x=np.zeros(dx.shape[0])
y=np.zeros(dx.shape[0])

norm=max(dx)
theta=cour*math.pi/2.
for i in range(err.shape[0]) :
  err[i]=math.log(err[i],2)
for i in range(dx.shape[0]) :
  r[i]=math.log(dx[i]/norm,2)+8
for i in range(theta.shape[0]) :
  x[i]=r[i]*math.cos(theta[i])
  y[i]=r[i]*math.sin(theta[i])

ngrid = 800*2
levels=[-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7] 

xi = np.linspace(0, 8, ngrid)
yi = np.linspace(0, 8, ngrid)
zi = griddata(x,y,err,xi,yi,interp='linear')

fig = plt.figure()
plt.contour(xi,yi,zi,levels,linewidths=0.5,colors='k')
plt.contourf(xi,yi,zi,levels,cmap=plt.cm.jet)
plt.colorbar() # draw colorbar

plt.xlim(0,8)
plt.ylim(0,8)
plt.title('log2(err)')
plt.xlabel('r; C=0')
plt.ylabel('r; C=1')

fig.savefig(sys.argv[1]+'.pdf')
fig.clf
