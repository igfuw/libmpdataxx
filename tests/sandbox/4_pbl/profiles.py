import h5py
import numpy as np
from sys import argv
import matplotlib.pyplot as plt

fname = argv[1]
h = h5py.File(fname, "r")

dz = 30. # TODO: get it from coords

w = h["w"]
tht = h["tht"]

_, _, NZ = w.shape
Z = np.array([dz * k for k in np.arange(NZ)])

av_tht = np.mean(tht, axis = (0, 1))
dtht = tht - av_tht

av_w = np.mean(w, axis = (0, 1))
dw = w - av_w

hflux = np.mean(dtht * dw, axis = (0,1))

zi = Z[np.argmin(hflux)]

Z = Z / zi

H0 = 0.01
g = 10.0
tht_ref = 300.
w_scale = (g / tht_ref * zi * H0) ** (1./3)
tht_scale = H0 / w_scale
t_scale = zi / w_scale

var_th = np.mean(dtht * dtht, axis = (0, 1))
var_w = np.mean(dw * dw, axis = (0, 1))

fig, axarr = plt.subplots(1, 3, figsize= (10, 6))

fig.subplots_adjust(top=0.85)
s = "$z_i = {zi}$ m, $w_* = {ws:.3}$ m/s, $t_* = {ts:.5}$ s, $T_* = {tht_s:.3}$ K"
fig.text(0.1, 0.95, s.format(zi = zi, ws = w_scale, ts = t_scale, tht_s = tht_scale), 
         bbox={'facecolor':'white', 'pad':5}
        )

axarr[0].plot(hflux / H0, Z, 'k-', lw = 1)
axarr[0].set_ylabel("$z / z_i$", fontsize = 20)
axarr[0].set_xlim([-0.5, 1.5])
axarr[0].set_ylim([0, 1.3])
axarr[0].set_title(r"$ \langle \Theta' w' \rangle / H_0$")

axarr[1].plot(var_th / (tht_scale ** 2), Z, 'k-', lw = 1)
axarr[1].set_xlim([0., 32.])
axarr[1].set_ylim([0, 1.3])
axarr[1].set_title(r"$ \langle \Theta' \Theta' \rangle / T_*^2$")

axarr[2].plot(var_w / (w_scale ** 2), Z, 'k-', lw = 1)
axarr[2].set_xlim([0., 0.5])
axarr[2].set_ylim([0, 1.3])
axarr[2].set_title(r"$ \langle w' w' \rangle / w_*^2$")

plt.savefig('profiles.pdf')
