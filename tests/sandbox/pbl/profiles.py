import h5py
import numpy as np
import os
from sys import argv
import matplotlib
matplotlib.use('Pdf')
import matplotlib.pyplot as plt

dirname = argv[1]
fnames = sorted(filter(lambda s : s.endswith('.h5'), os.listdir(dirname)))
fnames = [dirname + '/' + fn for fn in fnames]

# get timestep data
data_h = h5py.File(fnames[-1], "r")
w = data_h["w"]
tht = data_h["tht"]

# get const data
const_h = h5py.File(fnames[0], "r")
dz = const_h['advection'].attrs['dk'][0]
g = const_h['boussinesq'].attrs['g'][0]
tht_ref = const_h['boussinesq'].attrs['Tht_ref'][0] 
H0 = const_h['boussinesq'].attrs['hflux_const'][0] 

_, _, nz = w.shape
Z = np.array([dz * k for k in np.arange(nz)])

# calculate variances and the heat flux
av_tht = np.mean(tht, axis = (0, 1))
dtht = tht - av_tht
var_th = np.mean(dtht * dtht, axis = (0, 1))

av_w = np.mean(w, axis = (0, 1))
dw = w - av_w
var_w = np.mean(dw * dw, axis = (0, 1))

hflux = np.mean(dtht * dw, axis = (0,1))

# calculate the characteristic scales
zi = Z[np.argmin(hflux)]
w_scale = (g / tht_ref * zi * H0) ** (1./3)
tht_scale = H0 / w_scale
t_scale = zi / w_scale

# normalise profiles by the characteristic scales
Z /= zi
hflux /= H0
var_th /= tht_scale ** 2
var_w /= w_scale ** 2

# save profiles to a text file
outname = dirname[4:]
outfile = open('profiles_' + outname + '.txt', 'w')
for lev_d in zip(Z, hflux, var_th, var_w):
    outfile.write('{:6.2f} {:10.4e} {:10.4e} {:10.4e}\n'.format(*lev_d))

# plot profiles
fig, axarr = plt.subplots(1, 3, figsize= (10, 6))

fig.subplots_adjust(top=0.85)
s = "$z_i = {zi}$ m, $w_* = {ws:.3}$ m/s, $t_* = {ts:.5}$ s, $T_* = {tht_s:.3}$ K"
fig.text(0.1, 0.95, s.format(zi = zi, ws = w_scale, ts = t_scale, tht_s = tht_scale), 
         bbox={'facecolor':'white', 'pad':5}
        )

axarr[0].plot(hflux, Z, 'k-', lw = 1)
axarr[0].set_ylabel("$z / z_i$", fontsize = 20)
axarr[0].set_xlim([-0.5, 1.5])
axarr[0].set_ylim([0, 1.3])
axarr[0].set_title(r"$ \langle \Theta' w' \rangle / H_0$")

axarr[1].plot(var_th, Z, 'k-', lw = 1)
axarr[1].set_xlim([0., 32.])
axarr[1].set_ylim([0, 1.3])
axarr[1].set_title(r"$ \langle \Theta' \Theta' \rangle / T_*^2$")

axarr[2].plot(var_w, Z, 'k-', lw = 1)
axarr[2].set_xlim([0., 0.5])
axarr[2].set_ylim([0, 1.3])
axarr[2].set_title(r"$ \langle w' w' \rangle / w_*^2$")

plt.savefig('profiles_' + outname + '.pdf')
