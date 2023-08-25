"""

Explore the turbulent assumption in subglacial channels

"""

import numpy as np
import netCDF4 as nc

from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib import gridspec

import cmocean
from palettes.code import palettes,tools

import GladsPlot as gplt

# Set outputs to read in
out_fname = '../glads/01_kan_forcing/RUN/output_005_seasonal.nc'
nu = 1.79e-6
# mesh_fname = '../glads/data/mesh/mesh_04.nc'

# Read glads and mesh data
with nc.Dataset(out_fname, 'r') as out:
    # Channel fields
    Q = np.abs(out['Q'][:, :].data.T)
    S = np.abs(out['S_channel'][:, :].data.T)

    # Get floatation fraction.lege
    phi = out['phi'][:, :].data.T
    N = out['N'][:, :].data.T
    phi_0 = 9.81*1000*np.vstack(out['bed'][:].data)
    pw = phi - phi_0
    ff = pw/(N + pw)


v_channel = Q/S
radius_channel = np.sqrt(S/np.pi)
diam_channel = 2*radius_channel
Re_channel = diam_channel*v_channel/nu

# tindices = [365 + 125, 365 + 195]
# tindices = 365 + 125 + np.arange(0, 100, 1)
tindices = 365 + np.arange(0, 365)
# sort_index = np.argsort(Re_channel, axis=0)
# print(sort_index.shape)
# Re_sort = np.take_along_axis(Re_channel, sort_index, axis=0)
# Q_sort = np.take_along_axis(Q, sort_index, axis=0)
# Q_sort = Q[sort_index]
# Re_sort = Re_channel.flatten()[sort_index]
# Q_sort = Q.flatten()[sort_index]
# Q_cum = np.cumsum(Q)


bins = 10.**np.arange(-6, 8+1)
fig, ax = plt.subplots()
ax.hist(Re_channel.flatten(), bins=bins)

ax2 = ax.twinx()

for tindex in tindices:
    argsort = np.argsort(Re_channel[:, tindex])
    Re_sort = Re_channel[argsort, tindex]
    Q_sort = Q[argsort, tindex]
    ax2.plot(Re_sort, np.cumsum(Q_sort), color='k', alpha=0.1)
    
ax.set_xscale('log')
ax.set_yscale('log')

ax2.set_yscale('log')

ax.set_xlim([1e-6, 1e8])
ax.set_ylim([1e2, 1e8])

ax.axvline(2e3, color='k',linewidth=1)
ax.grid(linestyle='dashed',linewidth=0.5)

ax.set_xlabel(r"$\rm{Re}$")
ax.set_ylabel("Number of channels")

ax2.set_ylabel(r"Cumulative discharge $q~(\rm{m}^3~\rm{s}^{-1})$")

fig.subplots_adjust(right=0.85)

fig.savefig("figures/aux/channel_Re.png", dpi=600)

plt.show()


