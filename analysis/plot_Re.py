"""

Plot Reynolds number

"""

import numpy as np

import netCDF4 as nc
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.tri import Triangulation
import netCDF4 as nc

import cmocean
from palettes.code import palettes, tools

import GladsPlot as gplt

import helpers
import defaults


figsize=(7, 3.5)

gs_kwargs=dict(wspace=0.075, hspace=0.2, 
        left=0.09, right=0.95, bottom=0.12,
        top=0.95)

def plot_Re(fnames, figname, tslice=defaults.tslice,
    figsize=figsize, gs_kwargs=gs_kwargs, labels=defaults.labels,
    map_cmap=None, line_cmap=defaults.cmaps['Q'],
    Qmin=10, Qmax=100, Re_ylim=[0, 1e4], ff_yticks=None):

    ## CONFIG

    n_cases = len(fnames)

    if map_cmap is None:
        cm1 = palettes.get_cmap('BrownGray')
        cm2 = palettes.get_cmap('blue-8').reversed()
        z1 = int(200*2000/Re_ylim[1])
        z2 = 200 - z1
        print(z1)
        print(z2)
        map_cmap = tools.join_cmaps(cm1, cm2, N1=z1, N2=z2, average=10)

    ## Start the figure
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(3, 2, **gs_kwargs)
    cax_gs = gridspec.GridSpecFromSubplotSpec(4, 1,
        subplot_spec=gs[1, 1], height_ratios=(100, 30, 30, 100),
        hspace=0.2)

    # gridspec.GridSpecFromSubplotSpec(n_cases+2, 2, 
        # subplot_spec=global_gs[:, 0], width_ratios=(100, 4), height_ratios=hratios,
        # hspace=0.1  , wspace=0.1)

    ii_rows = [0, 0, 1, 2, 2]
    ii_cols = [0, 1, 0, 0, 1]
    axs = np.array([fig.add_subplot(gs[i,j]) for i,j in zip(ii_rows, ii_cols)])
    alphabet = ['a', 'b', 'c', 'd', 'e']


    # Start reading the data
    for ii in range(n_cases):
        fname = fnames[ii]
        print(fname)
        ax = axs[ii]

        out = nc.Dataset(fname)
        nodes = out['nodes'][:].data.T
        connect = out['connect'][:].data.T.astype(int) - 1
        connect_edge = out['connect_edge'][:].data.T.astype(int) - 1
        elements = out['elements'][:].data.T

        Q = np.abs(out['Q'][tslice, :].data.T)

        q = out['qs'][tslice, :].data.T
        qs = np.sqrt(q[:, 0]**2 + q[:, 1]**2)
        Re = qs/float(out['para/nu'][:].data)
        mtri = Triangulation(nodes[:, 0]/1e3, nodes[:, 1]/1e3, connect)

        Re_color = ax.tripcolor(mtri, Re[:],
            cmap=map_cmap, vmin=Re_ylim[0], vmax=Re_ylim[1])

        lc = gplt.plot_edge_data(nodes/1e3, connect_edge, Q[:],
            line_cmap, vmin=Qmin, vmax=Qmax)
        ax.add_collection(lc)


        ax.set_xlim([0, 100])
        ax.set_ylim([0, 25])
        ax.set_yticks([0, 12.5, 25])
        ax.set_aspect('equal')

        ax.text(0.95, 0.95, alphabet[ii], transform=ax.transAxes,
            va='top', ha='right', fontweight='bold')
        ax.text(0.95, 0.05, labels[ii], transform=ax.transAxes,
            va='bottom', ha='right', color='k', fontsize=8)
        if ii>=3:
            ax.set_xlabel('x (km)')
        else:
            ax.set_xticklabels([])

        if ii==1 or ii==4:
            ax.set_yticklabels([])
        
        if ii==2:
            ax.set_ylabel('y (km)')
        
        if ii==2:
            ax.plot([-1, 0], [-1, 0], color=(1, 1, 1, 0), label='  ')
    

    cax1 = fig.add_subplot(cax_gs[1])
    cax2 = fig.add_subplot(cax_gs[2])

    fig.colorbar(Re_color, cax=cax1, orientation='horizontal')
    cb2 = fig.colorbar(lc, cax=cax2, orientation='horizontal')

    cax1.set_yticks([])
    # cax1.tick_top()
    cax1.xaxis.tick_top()
    cax1.xaxis.set_label_position('top')
    cax1.set_xlabel(r'$\rm{Re}$')
    
    cax2.set_yticks([])
    cax2.set_xlabel(r'Q (m$^3$ s$^{-1}$)')

    cticks = np.linspace(0, Qmax, 6)
    cticks[0] = Qmin
    cticks = np.unique(cticks)
    cb2.set_ticks(cticks)
    fig.savefig(figname, dpi=600)
    return fig

# if __name__ == '__main__':
#     cases = [1, 2, 3, 4, 5]
#     fpattern = '../glads/_00_shmip_forcing_shmip_topo/RUN/output_%03d_seasonal.nc'
#     fnames = [fpattern%caseid for caseid in cases]
#     print(fnames)

#     figname = 'tmp.png'

#     plot_Re(fnames, figname, Qmin=10, Qmax=100, Re_ylim=(0, 4e3))
