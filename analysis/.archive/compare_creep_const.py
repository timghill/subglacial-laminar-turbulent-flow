"""
Compare steady results for different creep const values
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

import netCDF4 as nc

pattern_creep_glads = '../glads/00_shmip_forcing_shmip_topo/RUN/output_%03d_seasonal.nc'
cases_creep_glads = [201, 202, 203, 204, 205]
fnames_creep_glads = [pattern_creep_glads % caseid for caseid in cases_creep_glads]

pattern_creep_paterson = '../glads/00b_shmip_forcing_shmip_topo_creep_const/RUN/output_%03d_seasonal.nc'
cases_creep_paterson = [1, 2, 3, 4, 5]
fnames_creep_paterson = [pattern_creep_paterson % caseid for caseid in cases_creep_paterson]

xb = 30e3
band_width = 5e3

colors = np.array([[0.420, 0.510, 0.620, 1],
                   [0.579, 0.677, 0.781, 1],
                   [0.500, 0.500, 0.500, 1],
                   [0.859, 0.683, 0.275, 1],
                   [0.929, 0.835, 0.408, 1]])


fig = plt.figure(figsize=(7, 5))
gs = GridSpec(5, 2, height_ratios=(100, 100, 100, 100, 100))
axs = np.array([[fig.add_subplot(gs[i,j]) for j in range(2)] for i in range(5)])
# caxs = np.array([fig.add_subplot(gs[0, i]) for i in range(2)])
for ii in range(len(fnames_creep_glads)):
    ax1 = axs[ii, 0]
    ax2 = axs[ii, 1]

    out_glads = nc.Dataset(fnames_creep_glads[ii])
    out_paterson = nc.Dataset(fnames_creep_paterson[ii])

    N_glads = out_glads['N'][:].data.T
    N_paterson = out_paterson['N'][:].data.T

    nodes = out_glads['nodes'][:].data.T
    node_mask = np.logical_and(nodes[:, 0]>=(xb - band_width/2), nodes[:, 0]<=(xb + band_width/2))

    Q_glads = np.abs(out_glads['Q'][:].data.T)
    Q_paterson = np.abs(out_paterson['Q'][:].data.T)

    tt = out_glads['time'][:].data.T
    time = tt/365/86400 - 100

    ax1.plot(time, np.mean(N_glads[node_mask, :], axis=0)/1e6, color=colors[ii])
    ax2.plot(time, np.mean(N_paterson[node_mask, :], axis=0)/1e6, color=colors[ii])

    ax1.set_ylim([-3, 6])
    ax2.set_ylim([-3, 6])
    ax1.grid()
    ax2.grid()

    creep_glads = out_glads['para/creep_const_s'][:].data
    creep_paterson = out_paterson['para/creep_const_s'][:].data

    if ii==0:
        ax1.set_title(r'$\tilde A = %.2e$' % creep_glads)
        ax2.set_title(r'$\tilde A = %.2e$' % creep_paterson)

    ax1.set_ylabel('N (MPa)')

    if ii<4:
        ax1.set_xticklabels([])
        ax2.set_xticklabels([])
    
    N_mean_glads = np.mean(N_glads[node_mask, :], axis=0)
    N_mean_paterson = np.mean(N_paterson[node_mask, :], axis=0)
    ax1.text(0.95, 0.9, 'N = %.2f MPa' % (N_mean_glads[-1]/1e6), transform=ax1.transAxes, ha='right', va='top')
    ax2.text(0.95, 0.9, 'N = %.2f MPa' % (N_mean_paterson[-1]/1e6), transform=ax2.transAxes, ha='right', va='top')

ax1.set_xlabel('Year')
ax2.set_xlabel('Year')

fig.savefig('seasonal_creep_comparison.png', dpi=600)

plt.show()
