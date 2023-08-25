"""

Plot steady floatation fraction maps and Re

"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.tri import Triangulation
import netCDF4 as nc

import cmocean
from palettes.code import palettes, tools

import GladsPlot as gplt

import helpers
import defaults


figsize=(7, 6)

gs_kwargs=dict(wspace=0.05, hspace=0.2, 
        left=0.09, right=0.92, bottom=0.08,
        top=0.915, width_ratios=(100, 4, 100, 3))



def plot_pressure_maps_timeseries(fnames, figname, tslice=-1, 
    figsize=figsize, gs_kwargs=gs_kwargs, labels=defaults.labels, 
    colors=defaults.colors, map_cmap=defaults.cmaps['floatation'],
    line_cmap=defaults.cmaps['Q'], Qmin=10, Qmax=100,
    ff_ylim=[0, 1], ff_yticks=[0, 0.25, 0.5, 0.75, 1],
    Re_ylim=[0, 4000], lws=defaults.linewidths,
    linestyles=defaults.linestyles, zorders=defaults.zorders):
    """
    Plot 2D floatation fraction maps and timeseries.

    [More description]

    Inputs:
    --------
    fnames : Iterable of str
        List of paths to model outputs for a simulation case
    
    figname : str
        Path to save figure
    
    Options:
    --------
    tslice : int
        Time index for 2D pressure maps
    
    x_bands : Array-like of floats
        x distances in km for timeseries. Timeseries are computed for mean
        floatation fraction in bands [xb-band_width/2, xb+band_width/2]
        for xb in x_bands.Temp 
    
    band_width : float
        Width of bands for band-averaging
    
    figsize : (width, height) = (6, 7)
        Tuple in inches of figure size. 
    
    gs_Kwargs : dict
        Dictionary of options passed to gridspec.GridSpec for global config
    
    labels : list of str
        List of strings specifying labels for legend
    
    colors : (M, 3) or (M, 4) array-like of rgb or rgba values

    map_cmap : LinearSegmentedColormap 
        Colormap object for tripcolor panels
    
    line_cmap : LinearSegmentedColormap
        Colormap object for plotting channel discharge
    
    Qmin, Qmax : float
        Min and max discharge for channel discharge colorbar
    
    tlim : (tmin, tmax)
        x-axis bounds for timeseries panels
    
    t_ticks : Array-like1
        x-axis ticks for timeseries panels

    ff_ylim : (ymin, ymax) for ff panels
    
    Returns:
    --------
        fig : matplotlib figure object
    
    """ 
    ## CONFIG

    n_cases = len(fnames)

    ## Start the figure
    fig = plt.figure(figsize=figsize)

    # 2 columns: 5 maps with space for colorbars
    hratios = 100*np.ones(n_cases+2)
    hratios[0] = 10
    hratios[-1] = 150
    gs_maps = gridspec.GridSpec(n_cases+2, 4, height_ratios=hratios, **gs_kwargs)

    # Initialize axes
    axs_maps = np.array([[fig.add_subplot(gs_maps[i+1, 2*j]) for j in range(2)] for i in range(n_cases)])
    axs_scatter = np.array([fig.add_subplot(gs_maps[-1, 2*i]) for i in range(2)])

    # Set style for panel labels
    map_alphabet = ['a', 'b', 'c', 'd', 'e', 'f']
    text_args = {'fontweight':'bold'}

    # Start reading the data
    for ii in range(n_cases):
        fname = fnames[ii]
        print(fname)

        out = nc.Dataset(fname, 'r')
        
        walltime = float(out['model/wallclock'][:])
        print('Walltime (hours):', walltime/3600)

        nodes = out['nodes'][:].data.T
        elements = out['elements'][:].data.T
        connect = out['connect'][:].data.T.astype(int) - 1
        connect_edge = out['connect_edge'][:].data.T.astype(int) - 1

        # Channel fields
        Q = np.abs(out['Q'][:, :].data.T)

        # Get floatation fraction.lege
        phi = out['phi'][:, :].data.T
        N = out['N'][:, :].data.T
        phi_0 = 9.81*1000*np.vstack(out['bed'][:].data)
        pw = phi - phi_0
        ff = pw/(N + pw)

        qxy = out['qs'][:].data.T
        qs = np.sqrt(qxy[:, 0]**2 + qxy[:, 1]**2)
        nu = 1.79e-6
        omega = 1/2000
        Re = qs/nu

        tt = out['time'][:].data/86400/365 - 100

        # Initialize triangulation for faster plotting
        mtri = Triangulation(nodes[:, 0]/1e3, nodes[:, 1]/1e3, connect)
        
        # Map panel
        ax1 = axs_maps[ii, 0]
        ax2 = axs_maps[ii, 1]
        fcolor = ax1.tripcolor(mtri, ff[:, tslice], cmap=map_cmap, vmin=0, vmax=1)
        ax1.set_aspect('equal')
        ax1.set_xlim([0, 100])
        ax1.set_ylim([0, 25])

        ax2.set_xlim([0, 100])
        ax2.set_ylim([0, 25])

        
        cm1 = palettes.get_cmap('BrownGray')
        cm2 = palettes.get_cmap('blue-8').reversed()
        z1 = int(200*2000/Re_ylim[1])
        z2 = 200 - z1
        Re_map = tools.join_cmaps(cm1, cm2, N1=z1, N2=z2, average=10)

        Rcolor = ax2.tripcolor(mtri, Re[:, tslice], cmap=Re_map, vmin=Re_ylim[0], vmax=Re_ylim[1])

        # if ii==0:
        #     ax1.set_yticks([0, 12.5, 25])
        # else:
            # ax1.set_yticks([0, 12.5])
        ax1.set_yticks([0, 12.5, 25])

        lc = gplt.plot_edge_data(nodes/1e3, connect_edge, Q[:, tslice],
            line_cmap, vmin=Qmin, vmax=Qmax)
        lc2 = gplt.plot_edge_data(nodes/1e3, connect_edge, Q[:, tslice],
            line_cmap, vmin=Qmin, vmax=Qmax)
        ax1.add_collection(lc)
        ax2.add_collection(lc2)

        if ii<n_cases:
            ax1.set_xticklabels([])
            ax2.set_xticklabels([])
        
        ax2.set_yticklabels([])
        
        ax2.text(0.95, 0.95, map_alphabet[ii], transform=ax2.transAxes,
            va='top', ha='right', **text_args)
        ax2.text(0.95, 0.05, labels[ii], transform=ax2.transAxes,
            va='bottom', ha='right', fontsize=8, color='k')
        xmid, ff_avg = helpers.width_average(nodes, ff[:, tslice])

        xmid, Re_avg = helpers.width_average(elements, Re[:, tslice])

        ax_scatter = axs_scatter[0]
        ax_scatter.plot(xmid/1e3, ff_avg, color=colors[ii], label=labels[ii],
            linewidth=lws[ii], linestyle=linestyles[ii], zorder=zorders[ii])
        ax_scatter.set_ylim(ff_ylim)
        ax_scatter.set_yticks(ff_yticks)
        ax_scatter.set_xlim([0, 100])
        ax_scatter.grid(linestyle=':', linewidth=0.5)

        ax_scatter2 = axs_scatter[1]
        ax_scatter2.plot(xmid/1e3, Re_avg, color=colors[ii], label=labels[ii],
            linewidth=lws[ii], linestyle=linestyles[ii], zorder=zorders[ii])
        ax_scatter2.set_ylim(Re_ylim)
        # ax_scatter2.set_yticks(ff_yticks)
        ax_scatter2.set_xlim([0, 100])
        ax_scatter2.grid(linestyle=':', linewidth=0.5)

        if ii==0:
            ax_scatter.set_xlabel('x (km)')
            ax_scatter2.set_xlabel('x (km)')


        out.close()

    ax_scatter.set_ylabel(r'$p_{\rm{w}}/p_{\rm{i}}$')
    ax_scatter2.set_ylabel('Re')
    ax_scatter2.text(0.95, 0.95, map_alphabet[n_cases], transform=ax_scatter2.transAxes,
        va='top', ha='right', **text_args)
    ax_scatter2.yaxis.tick_right()
    ax_scatter2.yaxis.set_label_position('right')
    

    axs_maps[2, 0].set_ylabel('y (km)')

    cax1 = fig.add_subplot(gs_maps[0, 0])
    cax2 = fig.add_subplot(gs_maps[0, 2])
    cax3 = fig.add_subplot(gs_maps[1:6, 3])

    cb1 = fig.colorbar(fcolor, cax=cax1, orientation='horizontal')
    cb2 = fig.colorbar(Rcolor, cax=cax2, orientation='horizontal')
    cb3 = fig.colorbar(lc, cax=cax3)

    cax1.xaxis.tick_top()
    cax1.xaxis.set_label_position('top')

    cax2.xaxis.tick_top()
    cax2.xaxis.set_label_position('top')

    cticks = np.linspace(0, Qmax, 6)
    cticks[0] = Qmin
    cticks = np.unique(cticks)
    cax3.set_yticks(cticks)

    cb1.set_label(r'$p_{\rm{w}}/p_{\rm{i}}$')
    cb2.set_label(r'${\rm{Re}}$')
    cb3.set_label(r'$Q~(\rm{m}^3~\rm{s}^{-1})$', labelpad=0)

    fig.savefig(figname, dpi=600)
    return fig

if __name__=='__main__':
    ## Case 00: Flat topo, synthetic forcing
    cases = [1, 2, 3, 4, 5]
    fnames = ['../glads/00_synth_forcing/RUN/output_%03d_steady.nc'%caseid for caseid in cases]
    figname = 'figures/supplement/00_steady.png'
    fig_00 = plot_pressure_maps_timeseries(fnames, figname, Qmin=1, Qmax=100)

    ## Case 01: Alpha=4
    cases = [2, 3, 4, 5, 6]
    fnames = ['../glads/01_kan_forcing/RUN/output_%03d_steady.nc'%caseid for caseid in cases]
    figname = 'figures/aux/01_steady_alpha.png'
    cols = np.zeros((5, 4))
    cols[0:4] = defaults.colors[1:5]
    cols[-1] = [0, 0, 0, 1]

    lws = np.zeros(5)
    lws[:4] = defaults.linewidths[1:]
    lws[-1] = 1
    linestyles = ['solid', 'dashed', 'solid', 'solid', 'solid']
    zorders = [2, 5, 2, 2, 3]
    fig_01 = plot_pressure_maps_timeseries(fnames, figname, Qmin=1, Qmax=100,
        labels=['Turbulent 3/2', 'Laminar 3', 'Transition 5/4', 'Transition 3/2', 'Laminar 4'],
        lws=lws, linestyles=linestyles, zorders=zorders, colors=cols)

    ## Case S01b: Flat topo, synthetic forcing
    cases = [1, 2, 3, 4, 5]
    fnames = ['../glads/S01b_parameter_sensitivity/RUN/output_%03d_steady.nc'%caseid for caseid in cases]
    figname = 'figures/aux/S01b_steady.png'
    fig_00 = plot_pressure_maps_timeseries(fnames, figname, Qmin=1, Qmax=100)


    plt.show()
