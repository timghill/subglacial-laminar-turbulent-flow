"""

Plot floatation fraction maps and timeseries.

Function plot_pressure_maps_timeseries does the work,
this can be called from external scripts

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
    width_ratios = (100, 30, 110), 
    left=0.09, right=0.98, bottom=0.08,
    top=0.915)



def plot_pressure_maps_timeseries(fnames, figname, tslice=defaults.tslice, 
    x_bands=defaults.x_bands, band_width=defaults.band_width, 
    figsize=figsize, gs_kwargs=gs_kwargs, labels=defaults.labels, 
    colors=defaults.colors, map_cmap=defaults.cmaps['floatation'],
    line_cmap=defaults.cmaps['Q'], Qmin=10, Qmax=100,
    t_lim=[1, 2], t_ticks=[1.0, 1.25, 1.5, 1.75, 2], ff_ylim=[0, 1.5],
    t_ticklabels=None, t_xlabel='Year', ff_yticks=[0, 0.5, 1, 1.5],
    melt_forcing='SHMIP', fill_between=False,
    lws=defaults.linewidths, linestyles=defaults.linestyles,
    zorders=defaults.zorders):
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

    # Sort out melt forcing
    if melt_forcing=='KAN':
        tt_temp = np.loadtxt('../glads/data/kan_l_melt/KAN_L_2014_temp_clipped.txt', delimiter=',')
        tt_days = tt_temp[:, 0]
        temp_sl = tt_temp[:, 1]
        lr = -0.005
        DT = lr*390
        temp_fun = lambda t: np.maximum(0*t, DT + np.interp(t%1, tt_days/365, temp_sl, left=0, right=0))
    elif melt_forcing=='KANadj':
        tt_temp = np.loadtxt('../glads/data/kan_l_melt/KAN_L_2014_temp_adjusted.txt', delimiter=',')
        tt_days = tt_temp[:, 0]
        temp_sl = tt_temp[:, 1]
        lr = -0.0075
        DT = lr*390
        temp_fun = lambda t: np.maximum(0*t, DT + np.interp(t%1, tt_days/365, temp_sl, left=0, right=0))
    elif melt_forcing=='SHMIP':
        temp_fun = lambda t: np.maximum(-16*np.cos(2*np.pi*t) - 5, 0*t)
    elif melt_forcing=='SHMIPadj':
        a = 9.0684
        DT = 390*0.0075
        DT_term = 0.005*390
        temp_fun = lambda t: -a*np.cos(2*np.pi*t) + a*(DT - 5)/16 - DT_term


    ## Start the figure
    fig = plt.figure(figsize=figsize)

    # A global gridspec giving two columns to work with
    global_gs = gridspec.GridSpec(1, 3, **gs_kwargs)

    # Left column: 3 timeseries panels and melt forcing
    hratios = 100*np.ones(len(x_bands)+2)
    hratios[0] = 8
    gs_timeseries = gridspec.GridSpecFromSubplotSpec(len(x_bands) + 2, 1, 
        subplot_spec=global_gs[:, 2], height_ratios=hratios)

    # Right column: 5 maps with space for colorbars
    hratios = 100*np.ones(n_cases+2)
    hratios[0] = 8
    hratios[-1] = 150
    gs_maps = gridspec.GridSpecFromSubplotSpec(n_cases+2, 2, 
        subplot_spec=global_gs[:, 0], width_ratios=(100, 4), height_ratios=hratios,
        hspace=0.1  , wspace=0.1)

    # Initialize axes
    axs_timeseries = np.array([fig.add_subplot(gs_timeseries[i+1, 0]) for i in range(len(x_bands) + 1)])
    axs_maps = np.array([fig.add_subplot(gs_maps[i+1, 0]) for i in range(n_cases)])
    ax_scatter = fig.add_subplot(gs_maps[-1, 0])

    # Set style for panel labels
    time_alphabet = ['g', 'h', 'i', 'j']
    map_alphabet = ['a', 'b', 'c', 'd', 'e', 'f']
    text_args = {'fontweight':'bold'}
    # Start reading the data
    for ii in range(n_cases):
        fname = fnames[ii]
        print(fname)

        with nc.Dataset(fname, 'r') as out:
            walltime = float(out['model/wallclock'][:])
            print('Walltime (hours):', walltime/3600)

            nodes = out['nodes'][:].data.T
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

            tt = out['time'][:].data/86400/365 - 100

        with nc.Dataset('../glads/data/mesh/mesh_04.nc', 'r') as dmesh:
            node_area = dmesh['tri/area_nodes'][:].data

        # Initialize triangulation for faster plotting
        mtri = Triangulation(nodes[:, 0]/1e3, nodes[:, 1]/1e3, connect)
        
        # Map panel
        mapax = axs_maps[ii]
        fcolor = mapax.tripcolor(mtri, ff[:, tslice], cmap=map_cmap, vmin=0, vmax=1)
        mapax.set_aspect('equal')
        mapax.set_xlim([0, 100])
        mapax.set_ylim([0, 25])

        # if ii==0:
        #     mapax.set_yticks([0, 12.5, 25])
        # else:
            # mapax.set_yticks([0, 12.5])
        mapax.set_yticks([0, 12.5, 25])

        lc = gplt.plot_edge_data(nodes/1e3, connect_edge, Q[:, tslice],
            line_cmap, vmin=Qmin, vmax=Qmax)
        mapax.add_collection(lc)

        if ii<n_cases:
            mapax.set_xticklabels([])
        
        mapax.text(0.95, 0.95, map_alphabet[ii], transform=mapax.transAxes,
            va='top', ha='right', **text_args)
        mapax.text(0.95, 0.05, labels[ii], transform=mapax.transAxes,
            va='bottom', ha='right', fontsize=8, color='w')
        xmid, ff_avg = helpers.weighted_width_average(nodes, ff[:, tslice], node_area)

        quantile_95 = lambda x: np.quantile(x, 0.95)

        if fill_between:
            ax_scatter.fill_between(xmid/1e3, ff_lower, ff_upper, facecolor=colors[ii], alpha=0.33,
                edgecolor=None, linewidth=linewidths[ii], zorder=zorders[ii],
                linestyle=linestyles[ii])
        ax_scatter.plot(xmid/1e3, ff_avg, color=colors[ii], label=labels[ii],
            linewidth=lws[ii], linestyle=linestyles[ii],
            zorder=zorders[ii])
        ax_scatter.set_ylim(ff_ylim)
        ax_scatter.set_yticks(ff_yticks)
        ax_scatter.set_xlim([0, 100])
        ax_scatter.grid(linestyle=':', linewidth=0.5)

        if ii==0:
            ax_scatter.set_xlabel('x (km)')

        # Timeseries
        for j, xb in enumerate(x_bands):
            xmin = xb - band_width/2
            xmax = xb + band_width/2
            node_mask = np.logical_and(nodes[:, 0]/1e3>=xmin, nodes[:, 0]/1e3<xmax)
            f_mean = np.sum(ff[node_mask, :]*np.vstack(node_area[node_mask]), axis=0)/np.sum(node_area[node_mask])
            f_lower = np.quantile(ff[node_mask, :], 0.025, axis=0)
            f_upper = np.quantile(ff[node_mask, :], 0.975, axis=0)
            timeax = axs_timeseries[j]

            if fill_between:
                timeax.fill_between(tt, f_lower, f_upper, facecolor=colors[ii], alpha=0.3)

            timeax.plot(tt, f_mean, label=labels[ii],
                color=colors[ii], linewidth=lws[ii], zorder=zorders[ii],
                linestyle=linestyles[ii])

            mapax.axvline(xb, color='w', linewidth=0.5)
            timeax.axvline(tslice/365, color='k', linewidth=0.5)

            timeax.text(0.025, 0.95, time_alphabet[j], transform=timeax.transAxes,
            va='top', ha='left', **text_args)

            timeax.set_xticklabels([])


    ax_scatter.set_ylabel(r'$p_{\rm{w}}/p_{\rm{i}}$')
    ax_scatter.text(0.95, 0.95, map_alphabet[n_cases], transform=ax_scatter.transAxes,
        va='top', ha='right', **text_args)


    melt = temp_fun(tt)
    melt_ax = axs_timeseries[-1]
    melt_ax.plot(tt, melt, color='k', linewidth=1)

    melt_ax.set_ylabel('Temperature ($^\circ$C)')
    melt_ax.set_ylim([0, 12])
    melt_ax.set_yticks([0, 4, 8, 12])
    melt_ax.grid(linestyle=':', linewidth=0.5)
    melt_ax.axvline(tslice/365, color='k', linewidth=0.5)
    melt_ax.text(0.025, 0.95, time_alphabet[j+1], transform=melt_ax.transAxes,
        va='top', ha='left', **text_args)

    axs_maps[2].set_ylabel('y (km)')

    cax1 = fig.add_subplot(gs_maps[1:6, 1])
    cax2 = fig.add_subplot(gs_maps[0, 0])

    fig.colorbar(fcolor, cax=cax1, extend='max')
    fig.colorbar(lc, cax=cax2, orientation='horizontal', extend='both')

    cax1.xaxis.tick_top()
    cax1.xaxis.set_label_position('top')

    cticks = np.linspace(0, Qmax, 6)
    cticks[0] = Qmin
    cticks = np.unique(cticks)
    print(Qmin)
    print(Qmax)
    print(cticks)
    cax2.set_xticks(cticks)

    cax2.xaxis.tick_top()
    cax2.xaxis.set_label_position('top')

    cax1.text(0, 1.1, r'$p_{\rm{w}}/p_{\rm{i}}$', transform=cax1.transAxes)
    cax2.set_xlabel(r'$Q~(\rm{m}^3~\rm{s}^{-1})$')

    axs_timeseries[0].legend(bbox_to_anchor=[0, 1.02, 1., 0.102], loc='lower left',
        ncol=2, mode='expand', borderaxespad=0.05, frameon=False, borderpad=0)
    for j, xb in enumerate(x_bands):
        axi = axs_timeseries[j]
        axi.set_xticklabels([])
        axi.set_xlim(t_lim)
        axi.set_ylim(ff_ylim)
        axi.set_yticks(ff_yticks)
        axi.set_xticks(t_ticks)
        axi.grid(linestyle=':', linewidth=0.5)

        ax_scatter.axvline(xb, color='k', linewidth=0.5)

    melt_ax.set_xlim(t_lim)
    melt_ax.set_xticks(t_ticks)
    if t_ticklabels:
        melt_ax.set_xticklabels(t_ticklabels)

    axs_timeseries[-1].set_xlabel(t_xlabel)
    axs_timeseries[1].set_ylabel(r'$p_{\rm{w}}/p_{\rm{i}}$')

    fig.savefig(figname, dpi=600)
    return fig
