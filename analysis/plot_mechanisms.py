"""

Investigate mechanisms behind differences for turbulent/laminar models

"""

import numpy as np
import netCDF4 as nc

from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

from scipy import interpolate

import defaults
from helpers import weighted_width_average

figsize=(6, 6)

def plot_mechanisms(fnames, figname, models, tslice=defaults.tslice, 
    x_band=defaults.x_bands[1], band_width=defaults.band_width, 
    figsize=figsize, labels=defaults.labels, 
    colors=defaults.colors, tlim=[3, 9], t_ticks=[4, 5, 6, 7, 8, 9],
    t_ticklabels=['May', '', 'July', '', 'Sep', ''],
    k_turb=1, rhow=1000, g=9.81, nu=1.79e-6, omega=1/1000,
    lws=defaults.linewidths, linestyles=defaults.linestyles,
    zorders=defaults.zorders):
    """Plot components of discharge parameterization to attribute differences.

    Inputs:
    -------
    fnames : list of str
        List of filenames pointing to model outputs
    
    figname : str
        Path to save figure
    
    models : list of str
        For each fname, specify which parameterization was used. One of:
            'Turbulent 5/4'
            'Turbulent 3/2'
            'Laminar'
            'Transition 5/4'
            'Transition 3/2'
    
    x_band : float
        center of band for spatial averaging in km
    
    band_width : float
        width of band for spatial averaging in km
    
    figsize : (width, height) in inches
    
    labels : list of str
            list of labels for legend
    
    colors : (M, 3) or (M, 4) array of rgb or rgba

    tlim : (tmin, tmax)
        x-axis bounds for timeseries panels
    
    t_ticks : array-like of floats
        x-axis tick locations for timeseries panels
    
    k_turb : float
        Turbulent conductivity for k_eff/k_turb panels
    
    rhow : float density of water

    g : float gravity

    nu : float kinematic viscosity of water at 0C

    omega : transition parameter

    Returns:
    -------
    fig : figure object
    """

    textx = 0.05
    texty = 0.95
    textfmt = {'fontweight':'bold', 'fontsize':10,
                'ha':'left', 'va':'top'}

    # Get filenames and num cases4
    n_cases = len(fnames)

    tband_min = x_band - band_width/2
    tband_max = x_band + band_width/2

    # Set figure
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(5, 2, hspace=0.15, top=0.9, bottom=0.075, 
        right=0.975, left=0.125, wspace=0.1)
    alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
    axs = np.array([[fig.add_subplot(gs[ii, jj]) for jj in range(2)] for ii in range(5)])

    labels = ['Turbulent 5/4', 'Turbulent 3/2', 'Laminar', ' ', 'Transition 5/4', 'Transition 3/2']

    w = 25e3    # Domain width (km)
    for ii in range(n_cases):
        fname = fnames[ii]
        print(fname)

        out = nc.Dataset(fname, 'r')

        nodes = out['nodes'][:].data.T
        elements = out['elements'][:].data.T
        connect_edge = out['connect_edge'][:2, :].data.T.astype(int) - 1

        band_mask = np.logical_and(elements[:, 0]/1e3>=tband_min, elements[:, 0]/1e3<=tband_max)

        n_nodes = nodes.shape[0]
        n_elements = elements.shape[0]
        n_edges = connect_edge.shape[0]

        time = out['time'][:].data.T
        tt = 12*(time/86400/365 - 101)
        n_time = len(time)

        ## Net sheet and channel flux
        qs_xy = out['qs'][:].data.T
        qs = np.sqrt(qs_xy[:, 0]**2 + qs_xy[:, 1]**2)

        with nc.Dataset('../glads/data/mesh/mesh_04.nc', 'r') as dmesh:
            node_area = (dmesh['tri/area_nodes'][:].data)
            el_area = (dmesh['tri/area'][:].data)

        ## 1  Turbulence index (omega * Re)
        ax1 = axs[0, 0]
        ax2 = axs[0, 1]

        ## omega * Re
        lineargs = dict(linewidth=lws[ii], linestyle=linestyles[ii],
                        zorder=zorders[ii], color=colors[ii])
        xgrid, q_mean = weighted_width_average(elements, qs[:, tslice], el_area)
        Re_mean = q_mean/nu
        ax1.plot(xgrid/1e3, omega*Re_mean, **lineargs)
        ax1.grid()
        ax1.set_xlim([0, 100])
        # ax1.set_ylim([0, 12])
        ax1.set_ylim([1e-4, 50])
        # ax1.set_yticks([0, 1, 5, 10])
        ax1.set_ylabel(r'$\omega{\rm{Re}}$', labelpad=4)
        # ax1.legend(labels=labels, bbox_to_anchor=[0., 1, 2.2, 0.2], ncol=3,
        #     loc='lower center', frameon=False, mode='expand')
        ax1.text(textx, texty, 'a', transform=ax1.transAxes, **textfmt)
        ax1.set_yscale('log')

        if ii==2:
            ax1.plot([-1, 0], [-1, 0], color=(1, 1, 1, 0))

        Re = qs/nu
        mean_ts = np.average(Re[band_mask], axis=0, weights=el_area[band_mask])
        ax2.plot(tt, omega * mean_ts, **lineargs)
        ax2.grid()
        ax2.set_ylim([1e-4, 50])
        # ax2.set_yticks([0, 1, 5, 10])
        ax2.text(textx, texty, 'b', transform=ax2.transAxes, **textfmt)
        ax2.set_yscale('log')
        
        ## Compute gradphi for each case
        
        hs = out['h_sheet'][:].data.T
        interpolator = interpolate.LinearNDInterpolator(nodes, hs)
        hs_element = interpolator(elements)
        k = float(out['para/cond_s'][:].data)

        if models[ii]=='Turbulent 5/4':
            gradphi = (qs/k/hs_element**(5/4))**2
        elif models[ii]=='Turbulent 3/2':
            gradphi = (qs/k/hs_element**(3/2))**2
        elif models[ii]=='Laminar':
            gradphi = (qs/k/hs_element**(3))**1
        elif models[ii]=='Laminar 4':
            gradphi = (qs/k/hs_element**4)**1
        elif models[ii]=='Transition 5/4':
            h_bed = float(out['para/h_bed'][:].data)
            gradphi = ((qs + omega/nu * (hs_element/h_bed)**(1/2) * qs**2)/k/hs_element**(3))**1
        elif models[ii]=='Transition 3/2':
            gradphi = ((qs + omega/nu * qs**2)/k/hs_element**(3))**1

        T = rhow*g*(qs/gradphi)
        
        ## Transmissivity
        x_mean, T_mean = weighted_width_average(elements, T[:, tslice], el_area)
        ax1 = axs[1, 0]
        ax1.plot(x_mean/1e3, T_mean, **lineargs)
        # ax1.set_ylim([0, 2.5])
        ax1.set_ylabel(r'T (m$^2$ s$^{-1}$)', labelpad=4)
        ax1.text(textx, texty, 'c', transform=ax1.transAxes, **textfmt)
        ax1.set_yscale('log')
        ax1.set_ylim([1e-4, 1e1])
        ax1.set_yticks([1e-4, 1e-2, 1e0])

        T[elements[:, 0]<10e3, :] = np.nan
        ax2 = axs[1, 1]
        ax2.plot(tt, np.average(T[band_mask, :], axis=0, weights=el_area[band_mask]), **lineargs)
        ax2.text(textx, texty, 'd', transform=ax2.transAxes, **textfmt)
        ax2.set_yscale('log')
        ax2.set_ylim([1e-4, 1e1])
        ax2.set_yticks([1e-4, 1e-2, 1e0])

        ## Water thickness
        ax1 = axs[2, 0]
        x_mean, hs_mean = weighted_width_average(elements, hs_element[:, tslice], el_area)
        ax1.plot(x_mean/1e3, hs_mean, **lineargs)
        ax1.set_ylabel(r'$h$ (m)', labelpad=4)
        ax1.text(textx, texty, 'e', transform=ax1.transAxes, **textfmt)
        ax1.set_yscale('log')
        ax1.set_ylim([1e-4, 1e1])
        
        ax2 = axs[2, 1]
        ax2.plot(tt, np.average(hs_element[band_mask, :], axis=0, weights=el_area[band_mask]), **lineargs)
        ax2.set_yscale('log')
        ax2.text(textx, texty, 'f', transform=ax2.transAxes, **textfmt)
        ax2.set_ylim([1e-4, 1e1])

        # Potential gradient
        ax1 = axs[3, 0]
        x_mean, gradphi_mean = weighted_width_average(elements, gradphi[:, tslice], el_area, dx=2e3)
        ax1.plot(x_mean/1e3, gradphi_mean, **lineargs)
        ax1.set_ylabel(r'$|\nabla \phi|$ (Pa m$^{-1}$)')
        ax1.text(textx, texty, 'g', transform=ax1.transAxes, **textfmt)
        ax1.set_ylim([0, 500])

        ax2 = axs[3, 1]
        ax2.plot(tt, np.average(gradphi[band_mask, :], axis=0, weights=el_area[band_mask]), **lineargs)
        ax2.set_ylim([0, 500])
        ax2.text(textx, texty, 'h', transform=ax2.transAxes, **textfmt)

        ## Effective turbulent conductivity
        k_eff = qs/hs_element**(5/4)/gradphi**(1/2)

        print('Effective turbulent conductivity:')

        print('Total variation:', k_eff.max()/k_eff.min())
    
        ax1 = axs[4, 0]
        x_mean, k_eff_mean = weighted_width_average(elements, k_eff[:, tslice], el_area)

        print('Spatial variation:', k_eff_mean.max()/k_eff_mean.min())
        ax1.plot(x_mean/1e3, k_eff_mean/k_turb, **lineargs)
        ax1.set_ylabel(r'$k_{\rm{eff}}/k_{\rm{turb}}$', labelpad=4)
        ax1.text(textx, texty, 'i', transform=ax1.transAxes, **textfmt)
        ax1.set_yscale('log')
        ax1.set_ylim([1e-2, 1e1])
        ax1.set_yticks([1e-2, 1e-1, 1e0, 1e1])

        k_eff_timeseries = np.average(k_eff[band_mask, :], axis=0, weights=el_area[band_mask])/k_turb
        print('Temporal variation:', k_eff_timeseries.max()/k_eff_timeseries.min())
        ax2 = axs[4, 1]
        ax2.plot(tt, k_eff_timeseries, **lineargs)
        ax2.grid()
        ax2.text(textx, texty, 'j', transform=ax2.transAxes, **textfmt)
        ax2.set_yscale('log')
        ax2.set_ylim([1e-2, 1e1])
        ax2.set_yticks([1e-2, 1e-1, 1e0, 1e1])


    ax1 = axs[0, 0]
    # ax1.plot([-1, 0], [-1, 0], color='k')
    ax1.legend(labels=labels, bbox_to_anchor=[0., 1, 2, 0.2], ncol=3,
        loc='lower center', frameon=False, mode='expand')

    
    axs[-1, 0].set_xlabel('x (km)')
    axs[-1, 1].set_xlabel('Month')

    for ax in axs[:, 0]:
        ax.grid(True, linewidth=0.5, linestyle=':')

        ax.set_xlim([0, 100])
        ax.set_xticks([0, 20, 40, 60, 80, 100])
        ax.set_xticklabels([])
        ax.axvline(30, color='k', linewidth=0.5)

    ax.set_xticklabels(['0', '20', '40', '60', '80', '100'])

    for ax in axs[:, 1]:
        ax.grid(True, linewidth=0.5, linestyle=':')
        ax.set_xlim(tlim)
        ax.set_xticks(t_ticks)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.axvline((tslice/365 - 1)*12, color='k', linewidth=0.5)
    # ax.set_xticklabels([str(tii) for tii in t_ticks])
    ax.set_xticklabels(t_ticklabels)

    fig.savefig(figname, dpi=600)
    return fig

        
