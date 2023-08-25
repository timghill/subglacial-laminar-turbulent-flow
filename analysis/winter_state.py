"""

Explore the periodic winter state

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
mesh_fname = '../glads/data/mesh/mesh_04.nc'
tslice = 365 + 125
Qmin = -5
Qmax = 1

Smin = -5
Smax = 1

Q_cmap = palettes.get_cmap('BrownYellow')
S_cmap = palettes.get_cmap('BrownYellow')

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

with nc.Dataset(mesh_fname, 'r') as dmesh:
    nodes = dmesh['tri/nodes'][:].data.T
    connect = dmesh['tri/connect'][:].data.T.astype(int) - 1
    connect_edge = dmesh['tri/connect_edge'][:].data.T.astype(int) - 1

print('Max winter channel discharge:', np.max(Q[:, tslice]))
print('Max winter channel area:', np.max(S[:, tslice]))
Swinter = S[:, tslice]
n_channels = len(Swinter[Swinter>1e-3])
print('Number of channels with S>1e-3:', n_channels)

fig = plt.figure(figsize=(6, 4))
gs = gridspec.GridSpec(3, 2, width_ratios=(100, 3), 
    height_ratios=(10, 100, 100),
    hspace=0.25, wspace=0.05, left=0.12, right=0.9,
    bottom=0.15, top=0.85)
axs = np.array([fig.add_subplot(gs[i+1,0]) for i in range(2)])
cax1 = fig.add_subplot(gs[0,0])
cax2 = fig.add_subplot(gs[1, 1])
cax3 = fig.add_subplot(gs[2, 1])

mtri = Triangulation(nodes[:, 0]/1e3, nodes[:, 1]/1e3, connect)
pc = axs[0].tripcolor(mtri, ff[:, tslice], vmin=0, vmax=1, cmap=cmocean.cm.dense)
lc_Q = gplt.plot_edge_data(nodes/1e3, connect_edge, np.log10(Q[:, tslice]+1e-16),
    Q_cmap, vmin=Qmin, vmax=Qmax)
axs[0].add_collection(lc_Q)

_ = axs[1].tripcolor(mtri, ff[:, tslice], vmin=0, vmax=1, cmap=cmocean.cm.dense)
lc_S = gplt.plot_edge_data(nodes/1e3, connect_edge, np.log10(S[:, tslice]+1e-16),
    S_cmap, vmin=Smin, vmax=Smax)
axs[1].add_collection(lc_S)

cbf = fig.colorbar(pc, cax=cax1, orientation='horizontal')
cbQ = fig.colorbar(lc_Q, cax=cax2)
cbS = fig.colorbar(lc_S, cax=cax3)

cbf.set_label(r'$p_{\rm{w}}/p_{\rm{i}}$')
cbQ.set_label(r'$\log_{10} Q~(\rm{m}^3~\rm{s}^{-1})$')
cbS.set_label(r'$\log_{10} S~(\rm{m}^2)$')

cbQ.set_ticks([-5, -4, -3, -2, -1, 0, 1])
cbS.set_ticks([-5, -4, -3, -2, -1, 0, 1])

cax1.xaxis.tick_top()
cax1.xaxis.set_label_position('top')

alphabet = ['a', 'b']
for i,ax in enumerate(axs):
    ax.set_xlim([0, 100])
    ax.set_ylim([0, 25])
    ax.set_aspect('equal')
    ax.set_yticks([0, 12.5, 25])
    ax.text(0.95, 0.85, alphabet[i], fontweight='bold', color='k',
        transform=ax.transAxes)

axs[0].set_xticklabels([])
axs[-1].set_xlabel('Distance from terminus (km)')

fig.text(0.015, 0.5, 'Distance across glacier (km)',
    rotation='vertical', va='center')
fig.savefig('figures/aux/winter_channel_map.png', dpi=600)

hfig, hax = plt.subplots()
bins = (1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2)
hax.hist(Swinter, bins=bins)
hax.grid()
hax.set_xscale('log')
hax.set_yscale('log')
hax.set_xlabel('Channel area (m$^2$)')
hax.set_ylabel('Number of channels')

hfig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.15)
hfig.savefig('figures/aux/winter_channel_hist.png', dpi=600)

plt.show()
