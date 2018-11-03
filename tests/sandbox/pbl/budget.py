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
tht = np.array(data_h["tht"])
u = np.array(data_h["u"])
v = np.array(data_h["v"])
w = np.array(data_h["w"])
p = np.array(data_h["p"])

if "tke" in data_h:
    tke = np.array(data_h["tke"])
    subgrid = True
else:
    subgrid = False

# get const data
const_h = h5py.File(fnames[0], "r")
dt = const_h['advection'].attrs['dt'][0]
dz = const_h['advection'].attrs['dk'][0]
g = const_h['boussinesq'].attrs['g'][0]
tht_ref = const_h['boussinesq'].attrs['Tht_ref'][0]
H0 = const_h['boussinesq'].attrs['hflux_const'][0]
if subgrid:
    mix_len = np.array(const_h['boussinesq/mix_len'])

_, _, nz = w.shape
Z = np.array([dz * k for k in np.arange(nz)])

# calculate averages etc
av_tht = np.mean(tht, axis = (0, 1))
dtht = tht - av_tht

av_u = np.mean(u, axis = (0, 1))
du = u - av_u

av_v = np.mean(v, axis = (0, 1))
dv = v - av_v

av_w = np.mean(w, axis = (0, 1))
dw = w - av_w

if subgrid:
    av_tke = np.mean(tke, axis = (0, 1))

# heat flux
hflux = np.mean(dtht * dw, axis = (0,1))

# tke transport term
de = 0.5 * (du ** 2 + dv ** 2 + dw ** 2)
if subgrid:
    de += 2.0 / dt * (p - 2. / 3 * tke)
else:
    de += 2.0 / dt * p
tke_t_aux = np.mean(dw * de, axis = (0, 1))
tke_t = np.zeros(nz)
for k in range(1, nz - 1):
    tke_t[k] = (tke_t_aux[k + 1] - tke_t_aux[k - 1]) / (2 * dz)
tke_t[0] = (tke_t_aux[1] - tke_t_aux[0]) / dz
tke_t[nz - 1] = (tke_t_aux[nz - 1] - tke_t_aux[nz - 2]) / dz
tke_t *= -1

# dissipation of tke
if subgrid:
    ceps = 0.845
    dssp = -np.mean(ceps * tke * np.sqrt(tke) / mix_len, axis = (0, 1))
    dssp[0] = 0
else:
    dssp = np.zeros(nz)

# calculate the characteristic scales
zi = Z[np.argmin(hflux)]
w_scale = (g / tht_ref * zi * H0) ** (1./3)
tht_scale = H0 / w_scale
t_scale = zi / w_scale
Z = Z / zi
xnorm = zi / w_scale ** 3

# normalise profiles by the characteristic scales
hflux /= H0
tke_t *= xnorm
dssp *= xnorm

# budget inbalance
inb = dssp + hflux + tke_t

# save budget to a text file
outname = dirname[4:]
outfile = open('budget_' + outname + '.txt', 'w')
for lev_d in zip(Z, hflux, tke_t, dssp, inb):
    outfile.write('{:6.2f} {:10.1e} {:10.1e} {:10.1e} {:10.1e}\n'.format(*lev_d))

fig, axarr = plt.subplots(1, 1, figsize= (10, 6))
axarr.plot(hflux, Z, 'k-', lw = 1, label = 'B')
axarr.plot(tke_t, Z, 'k-.', lw = 1, label = 'T')
if subgrid:
    axarr.plot(dssp, Z, 'k--', lw = 1, label = 'D')
axarr.plot(inb, Z, 'k:', lw = 1, label = 'I')
axarr.set_ylabel("$z / z_i$", fontsize = 20)
axarr.set_xlim([-0.8, 1])
axarr.set_xticks([-0.8 + 0.2 * i for i in range(10)])
axarr.set_ylim([0, 1.3])
axarr.legend(fontsize = 20)

plt.savefig('budget_' + outname + '.pdf')
