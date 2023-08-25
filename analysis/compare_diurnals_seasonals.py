"""

Compare seasonal simulations with/without diurnals

"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
import netCDF4 as nc
import scipy.interpolate


cases = [1, 2, 3, 4]
no_diurnal_base = '../glads/01_kan_forcing/RUN/output_%03d_seasonal.nc'
yes_diurnal_base = '../glads/01c_kan_diurnal/RUN/output_%03d_seasonal.nc'
labels = ['Turbulent 5/4', 'Turbulent 3/2', 'Laminar',
    'Transition 5/4', 'Transition 3/2']

colors = np.array([[0.579, 0.677, 0.781, 1],
                   [0.199, 0.328, 0.492, 1],
                   [0.250, 0.250, 0.250, 1],
                   [0.929, 0.835, 0.408, 1],
                   [0.836, 0.590, 0.160, 1]])
linestyles = ['solid', 'dotted']
lws = [1.5, 1, 0.75, 1.5]
zorders = [1, 2, 4, 3]

fig, ax = plt.subplots(figsize=(7, 4))

for i,caseid in enumerate(cases):
    figname_no_diurnal = no_diurnal_base % caseid
    figname_yes_diurnal = yes_diurnal_base % caseid

    fnames = [figname_yes_diurnal, figname_no_diurnal]
    ffs = [[], []]
    tts = [[], []]
    figi, axi = plt.subplots()
    figj, axj = plt.subplots()
    for j in range(2):
        fname = fnames[j]
        print(fname)

        with nc.Dataset(fname, 'r') as out:
            phi = out['phi'][:].data.T
            N = out['N'][:].data.T
            bed = np.vstack(out['bed'][:].data)
            nodes = out['nodes'][:].data.T
            time = out['time'][:].data/86400/365 - 101
        
        phi_bed = 9.81*1000*bed
        pw = phi - phi_bed
        ff = pw/(N + pw)

        node_mask = np.logical_and(nodes[:, 0]<32.5e3, nodes[:, 0]>=27.5e3)
        ff_mean = np.mean(ff[node_mask, :], axis=0)
        axi.plot(time*12, ff_mean, color=colors[i], linestyle=linestyles[j])

        if linestyles[j]=='solid' and labels[i]=='Laminar':
            ls = 'dashed'
        else:
            ls = linestyles[j]
        if j==0:
            ax.plot(time*12, ff_mean, color=colors[i], label=labels[i],
                linewidth=lws[i], zorder=zorders[i], linestyle=ls)
        else:
            ax.plot(time*12, ff_mean, color=colors[i],
                linewidth=lws[i], zorder=zorders[i], linestyle=ls)
        ffs[j] = ff_mean
        tts[j] = time*12

    ff_no_diurnal = ffs[1]
    ff_yes_diurnal = ffs[0]
    ff_no_diurnal_interp_fun = scipy.interpolate.interp1d(tts[1], ffs[1], kind='linear')
    ff_no_diurnal_interp = ff_no_diurnal_interp_fun(tts[0])
    ff_diff = ff_yes_diurnal - ff_no_diurnal_interp
    axj.plot(tts[0], ff_diff, color=colors[i])
    axj.grid()
    axj.set_title(labels[i])
    axj.set_xlabel('Month')
    axj.set_ylabel(r'$\Delta$ flotation fraction')
    axj.set_xlim([5, 10])
    axj.set_ylim([-0.06, 0.06])

    axi.set_title(labels[i])
    axi.set_xlabel('Month')
    axi.set_ylabel(r'Floatation fraction ($p_{\rm{w}}/p_{\rm{i}}$)')
    axi.grid()
    # axi.set_ylim([0, 1])
    axi.set_xlim([5, 10])

    # figi.savefig('fig_diurnal_seasonal_%03d.png' % caseid, dpi=600)

    # figj.savefig('fig_diurnal_seasonal_diff_%03d.png' % caseid, dpi=600)

ax.grid()
ax.set_xlabel('Month')
ax.set_ylabel(r'Floatation fraction ($p_{\rm{w}}/p_{\rm{i}}$)')

ax.set_xlim([5, 6])
ax.legend()
fig.subplots_adjust(left=0.1, bottom=0.125, right=0.95, top=0.95)
fig.savefig('figures/supplement/00_compare_diurnal_seasonal.png', dpi=600)

# plt.show()
