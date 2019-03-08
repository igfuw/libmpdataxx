## @file
# @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
# @copyright University of Warsaw
# @date Januar 2012
# @section LICENSE
# GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)

# 1d isolines test from pks & wwg 1990

# for plotting in detached screen
import matplotlib
matplotlib.use('Svg')

from numpy import loadtxt, zeros, linspace
from sys import argv
from math import pi, log, cos, sin, sqrt
from matplotlib.mlab import griddata
from matplotlib.patches import Circle, Path, PathPatch

import matplotlib.pyplot as plt

dx, cour, err = loadtxt(argv[1], unpack=True)

mn=0
mx=8

theta = zeros(dx.shape[0])
r = zeros(dx.shape[0])
x = zeros(dx.shape[0])
y = zeros(dx.shape[0])

norm=max(dx)
theta = cour * pi / 2.
for i in range(err.shape[0]) :
  err[i] = log(err[i],2)
for i in range(dx.shape[0]) :
  r[i] = log(dx[i]/norm,2)+8
for i in range(theta.shape[0]) :
  x[i] = r[i] * cos(theta[i])
  y[i] = r[i] * sin(theta[i])

ngrid = 800 * 2
levels = range(-32,-3) #[-28,-27,-26,-25,-24,-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7]

xi = linspace(mn, mx, ngrid)
yi = linspace(mn, mx, ngrid)
zi = griddata(x,y,err,xi,yi,interp='linear')

fig = plt.gcf()
fig.gca().set_xlim(mn,mx)
fig.gca().set_ylim(mn,mx)
#fig.gca().set_title(r'log$_2$(err)', fontsize = 30)
fig.gca().set_xlabel('r; C=0', fontsize = 30)
fig.gca().set_ylabel('r; C=1', fontsize = 30)
fig.set_size_inches(9.78,8.2)
#fig.set_size_inches(8.2,8)
plt.tick_params(length=10, width=2, labelsize=24)


# colours
mpble = fig.gca().contourf(xi,yi,zi,levels,cmap=plt.cm.jet)
# colorbar
cbar = plt.colorbar(mpble)
#cbar.set_label('$log_2(err)$', labelpad=-40, padx = 10, y=0.45, fontsize=35)
cbar.set_label(r'log$_2$(err)', fontsize=30, labelpad = 20)
cbar.ax.tick_params(labelsize=24)
# grid
for r in range(mx,0,-1):
  zix = 1
  if (r==1): zix=10
  patch=PathPatch(
    Path(
      [(0,r),       (r*4*(sqrt(2)-1)/3,r), (r,r*4*(sqrt(2)-1)/3), (r,0)      ],
      [Path.MOVETO, Path.CURVE4,           Path.CURVE4,           Path.CURVE4]
    ),
    color='white',
    fill=(r==1),
    linewidth=1,
    zorder=zix
  )
  if (r!=1): fig.gca().add_patch(patch)
for i in range(0,19,3) :
  c = cos(theta[i])
  s = sin(theta[i])
  fig.gca().add_patch(
    PathPatch(
      Path(
	[(1*c,1*s),   (mx*c,mx*s)  ],
	[Path.MOVETO, Path.LINETO]
      ),
      color='white',
      fill=None,
      linewidth=1
    )
  )
# contours
fig.gca().contour(xi,yi,zi,levels,linewidths=1,colors='k')
# removing garbage near the origin
fig.gca().add_patch(patch)
# output
fig.savefig(argv[1]+'.svg')
