"""

Compare summer and winter water pressure statistics for turbulent, transition,
and laminar models

"""

import numpy as np
import netCDF4 as nc

from matplotlib import pyplot as plt

import defaults

# Global constants
rhow = 1000
rhoi = 910
g = 9.81

x_bands = np.array([15e3, 30e3, 70e3])
band_widths = 5e3
KAN_winter_tindex = [365 + 75, 365 + 125, 365 + 145, 365 + 145 + 150]
synth_winter_tindex = [365 + 50, 365 + 100, 365 + 120, 365 + 120 + 120]
labels = ['Turbulent 5/4', 'Turbulent 3/2',
          'Laminar',
          'Transition 5/4', 'Transition 3/2']
labels = [l.ljust(14) for l in labels]

def quantify_floatation(fnames, figname, tslices=None,
    x_bands=None, band_width=None, strfmt='%.4f', labels=None):
    """Docstring..."""
    if x_bands is None:
        ff_winters = np.zeros((len(fnames), 1))
        ff_summers = np.zeros((len(fnames), 1))
    else:
        ff_winters = np.zeros((len(fnames), len(x_bands)))
        ff_summers = np.zeros((len(fnames), len(x_bands)))
        ff_days = np.zeros((len(fnames), len(x_bands)))


    if tslices is None:
        tslices = [0, -1]

    if labels is None:
        labels = [fname.split('/')[-1] for fname in fnames]
    
    dmesh = nc.Dataset('../glads/data/mesh/mesh_04.nc')
    area_nodes = np.vstack(dmesh['tri/area_nodes'][:].data.T)
    fig, axs = plt.subplots(nrows=len(x_bands))
    # ff_strfmt = []
    # winter_strfmt = ''
    # summer_strfmt = ''
    for i in range(len(fnames)):
        # ff_strfmt.append(labels[i])
        out = nc.Dataset(fnames[i])
        phi = out['phi'][:].data.T
        N = out['N'][:].data.T
        
        bed = np.vstack(out['bed'][:].data)
        phi_elevation = rhow*g*bed

        pw = phi - phi_elevation
        ff = pw/(N + pw)
        nodex = np.vstack(out['nodes'][0].data)
        time = out['time'][:].data.T
        out.close()

        
        if x_bands is None:
            ff_global_mean = np.sum(ff*area_nodes)/np.sum(area_nodes)
            ff_winter = np.mean(ff_global_mean[tslices[0]:tslices[1]])
            ff_winters[i, 0] = ff_winter

            ff_summer = np.max(ff_global_mean)
            ff_summers[i, 0] = ff_summer

            # ff_strfmt.append(labels[i] + ':\t' + (strfmt % ff_mean))

        else:
            for k in range(len(x_bands)):
                xmin = x_bands[k] - band_width/2
                xmax = x_bands[k] + band_width/2
                mask = np.logical_and(nodex>=xmin, nodex<=xmax).flatten()

                ff_band_avg = np.sum(ff[mask]*(area_nodes[mask]), axis=0)/np.sum(area_nodes[mask])

                ff_winter = np.mean(ff_band_avg[tslices[0]:tslices[1]])
                # ff_summer = np.max(ff_band_avg)
                summer_slice = ff_band_avg[tslices[2]:tslices[3]]
                ff_summer = np.percentile(summer_slice, 95)

                days_overpressure = len(summer_slice[summer_slice>1])
                # print(ff_band_avg[tslices[2]:tslices[3]])
                ff_winters[i, k] = ff_winter
                ff_summers[i, k] = ff_summer
                ff_days[i, k] = days_overpressure
                
                tt = time/86400 - 100*365
                ax = axs[k]
                ax.plot(tt, ff_band_avg, color=defaults.colors[i])
                ax.axhline(ff_summer, linewidth=0.3, color=defaults.colors[i])
                ax.grid()
            # winter_strfmt = winter_strfmt + '\t'.join(ff_winters[i])
            # summer_strfmt = summer_strfmt + '\t'.join(ff_summers[i])
            # ff_strfmt.append(':\t' + winter_strfmt + '\n' + ''.ljust(14) + ':\t' + summer_strfmt)

        # ff_strfmt.append(labels[i] + ':\t' + (strfmt % ff_winter) + '\n' + ''.ljust(14) +  ':\t' + (strfmt % ff_summer))

    # ff_strfmt = '\n'.join(ff_strfmt)
    for k in range(len(x_bands)):
        ax = axs[k]
        ax.fill_betweenx([0, 2], [tslices[0], tslices[0]], [tslices[1], tslices[1]], color='gray', alpha=0.5)
        ax.fill_betweenx([0, 2], [tslices[2], tslices[2]], [tslices[3], tslices[3]], color='blue', alpha=0.2)
    fig.tight_layout()
    fig.savefig(figname, dpi=600)
    return (ff_winters, ff_summers, ff_days)

print('-------------------------------------------------')
print('    SYNTHETIC (case 00)')
print('-------------------------------------------------')
cases = [1, 2, 3, 4, 5]
synth_dir = '../glads/00_synth_forcing/RUN/output_%03d_seasonal.nc'
fnames = [synth_dir % caseid for caseid in cases]
figname = 'stats_synth.png'
winter, summer, days = quantify_floatation(fnames, figname, x_bands=x_bands, band_width=band_widths,
    tslices=synth_winter_tindex, labels=labels)
print('Winter:')
print(winter.T)
print('Summer:')
print(summer.T)
print('Days overpressure:')
print(days.T)


print('-------------------------------------------------')
print('    KAN (case 01)')
print('-------------------------------------------------')
cases = [1, 2, 3, 4, 5]
kan_dir = '../glads/01_kan_forcing/RUN/output_%03d_seasonal.nc'
fnames = [kan_dir % caseid for caseid in cases]
figname = 'stats_KAN.png'
winter, summer, days = quantify_floatation(fnames, figname, x_bands=x_bands, band_width=band_widths,
    tslices=KAN_winter_tindex, labels=labels)
print('Winter:')
print(winter.T)
print('Summer:')
print(summer.T)
print('Days overpressure:')
print(days.T)

print('-------------------------------------------------')
print('    SENSITIVITY')
print('-------------------------------------------------')
cases = [1, 2, 3, 4, 5]
kan_dir = '../glads/S01b_parameter_sensitivity/RUN/output_%03d_seasonal.nc'
fnames = [kan_dir % caseid for caseid in cases]
figname = 'stats_KAN_sensitivity.png'
winter, summer, days = quantify_floatation(fnames, figname, x_bands=x_bands, band_width=band_widths,
    tslices=KAN_winter_tindex, labels=labels)
print('Winter:')
print(winter.T)
print('Summer:')
print(summer.T)
print('Days overpressure:')
print(days.T)

plt.show()
