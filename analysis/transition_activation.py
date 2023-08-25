"""

Plot relative size of each term in the transition parameterization

"""

import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
from matplotlib import gridspec

import cmocean


def plot_activation(fname, figname, tlim=[3, 9]):
    fig = plt.figure(figsize=(6, 4))
    gs = gridspec.GridSpec(2, 2, width_ratios=(100, 50),
        height_ratios=(5, 100), hspace=0.05, wspace=0.15,
        bottom=0.125, top=0.875, right=0.975, left=0.1)
    
    ax0 = fig.add_subplot(gs[1, 0])
    ax1 = fig.add_subplot(gs[1, 1])
    cax = fig.add_subplot(gs[0, 0])

    dmesh = nc.Dataset('../glads/data/mesh/mesh_04.nc')
    nodes = dmesh['tri/nodes'][:].data.T
    connect = dmesh['tri/connect'][:].data.T
    # elements = dmesh['tri/elements'][:].data.T
    area = dmesh['tri/area'][:].data.T
    edges = dmesh['tri/edge_midpoints'][:].data.T

    with nc.Dataset(fname, 'r') as out:
        # h = out['h_sheet'][:].data.T
        qxy = out['qs'][:].data.T
        qs = np.sqrt(np.sum(qxy**2, axis=1))
        # omega = out['para/omega'][:].data
        omega = 1/2000
        k = out['para/cond_s'][:].data
        nu = out['para/nu'][:].data
        elements = out['elements'][:].data.T
        time = out['time'][:].data.T
        Q = np.abs(out['Q'][:].data.T)
    
    # Qlimit = np.max(edges[Q>1, 0], axis=0)
    # Qmask = Q>1
    # xvals = edges[:, 0][Qmask]
    xrep = np.zeros(Q.shape)
    print(xrep.shape)
    xrep[:, :] = np.vstack(edges[:, 0])
    xrep[Q<1] = np.nan
    Qlimit = np.nanmax(xrep, axis=0)
    xrep[Q<10] = np.nan
    Qlimit2 = np.nanmax(xrep, axis=0)

    
    term_laminar = qs
    term_turbulent = omega*(qs/nu)*qs

    # Width-average for plotting
    metric = lambda x: np.nanmean(x, axis=0)
    x, term_laminar_avg = width_average(elements, term_laminar, metric=metric)
    x, term_turbulent_avg = width_average(elements, term_turbulent, metric=metric)

    norm_laminar_avg = term_laminar_avg/(term_laminar_avg + term_turbulent_avg)
    norm_turbulent_avg = term_turbulent_avg/(term_laminar_avg + term_turbulent_avg)

    t_month = 12*(time - time[0])/365/86400 - 12
    [xx, tt] = np.meshgrid(x, t_month)
    pc = ax0.pcolormesh(xx/1e3, tt, norm_turbulent_avg.T, cmap=cmocean.cm.rain, vmin=0, vmax=1)
    Qlimit[np.isnan(Qlimit)] = 0
    ax0.plot(Qlimit/1e3, tt, color=(0.12, 0.15, 0.15), linestyle='solid', linewidth=0.2, alpha=1.)
    # ax0.plot(Qlimit2/1e3, tt, color=(0.12, 0.15, 0.15), linestyle='solid', linewidth=0.2, alpha=1.)

    fig.colorbar(pc, cax=cax, orientation='horizontal')

    norm_turbulent_domain_avg = np.mean(norm_turbulent_avg, axis=0)
    norm_turbulent_domain_max = np.max(norm_turbulent_avg, axis=0)
    ax1.plot(norm_turbulent_domain_avg, t_month, color='k', label='Mean')
    ax1.plot(norm_turbulent_domain_max, t_month, color='gray', label='Max')
    
    # Refine axes
    ax0.set_xlim([0, 100])
    ax1.set_xlim([0, 1])

    yticks = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    ax0.set_yticks(yticks)
    ax0.set_yticklabels(['Jan', 'Feb', 'Mar', 'April', 'May', 'June', 'July', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
    ax1.set_yticks(yticks)
    ax1.set_yticklabels([])
    ax0.set_ylim(tlim)
    ax1.set_ylim(tlim)

    ax0.set_xlabel('Upglacier distance (km)')
    ax0.set_ylabel('Month')
    ax1.set_xlabel(r'Fraction turbulent ($\gamma_{\rm{turb}}$)')
    ax1.grid()

    cax.xaxis.tick_top()
    cax.xaxis.set_label_position('top')
    cax.set_xlabel(r'Fraction turbulent ($\gamma_{\rm{turb}}$)')

    ax0.text(0.98, 0.95, 'a', transform=ax0.transAxes, va='top', ha='right', fontweight='bold')
    ax1.text(0.95, 0.95, 'b', transform=ax1.transAxes, va='top', ha='right', fontweight='bold')

    fig.savefig(figname, dpi=400)
    return fig


def width_average(xy, z, dx=1e3, metric=np.nanmean, xmin=0, xmax=100e3):
    """Width-average field z according to coordinates xy.

    Arguments:
    ---------
    xy : (N, 2) array in meters
    z : (N,) array
    dx : Flowline increment in meters
    area : Area corresponding to elements z is defined on
    metric : Function to evaluate (np.nanmean, np.median, etc...)

    Returns:
    --------
    xmid : Midpoint coordinates
    z_avg : z averaged over intervals centered on xmid

    """
    xedge = np.arange(xmin, xmax+dx, dx)
    xmid = 0.5*(xedge[1:] + xedge[:-1])

    xvec = np.reshape(xy[:, 0], (xy.shape[0], 1))
    gvec = np.array([xmid])

    tf_mask = np.logical_and(xvec>=(gvec-dx/2), xvec<=(gvec+dx/2))
    z_avg = np.zeros((len(xmid), z.shape[1]))
    for i in range(len(xmid)):
        z_avg[i] = metric(z[tf_mask[:, i]])
    
    return xmid, z_avg

if __name__=='__main__':
    fname = '../glads/00_synth_forcing/RUN/output_005_seasonal.nc'
    figname = 'figures/aux/00_activation.png'
    fig = plot_activation(fname, figname, tlim=[3, 9])

    fname = '../glads/01_kan_forcing/RUN/output_005_seasonal.nc'
    figname = 'figures/supplement/01_activation.png'
    fig = plot_activation(fname, figname, tlim=[4, 10])
    plt.show()
