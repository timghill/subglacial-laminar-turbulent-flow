"""

Call plot_Re to plot Reynolds number maps

"""

import numpy as np
from matplotlib import pyplot as plt

from plot_Re import plot_Re

## Case 00: Flat topo, synthetic forcing
cases = [1, 2, 3, 4, 5]
fnames = ['../glads/00_synth_forcing/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/main/00_Re_seasonal.png'
fig_00 = plot_Re(fnames, figname, Qmin=1, Qmax=100, Re_ylim=(0, 4e3), tslice=530)

## Case 00a: Flat topo, SHMIP forcing
cases = [1, 2, 3, 4, 5]
fnames = ['../glads/00a_shmip_forcing/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/00a_Re_seasonal_SHMIP_day515.png'
fig_00 = plot_Re(fnames, figname, Qmin=10, Qmax=200, Re_ylim=(0, 8e3), tslice=515)

## Case 00a: Flat topo, SHMIP forcing
cases = [1, 2, 3, 4, 5]
fnames = ['../glads/00a_shmip_forcing/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/00a_Re_seasonal_SHMIP_day530.png'
fig_00 = plot_Re(fnames, figname, Qmin=10, Qmax=200, Re_ylim=(0, 8e3), tslice=530)

## Case 01: Flat topo, KAN_L forcing
cases = [1, 2, 3, 4, 5]
pattern = '../glads/01_kan_forcing/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = 'figures/aux/01_Re_seasonal.png'
fig_01 = plot_Re(fnames, figname, Qmin=1, Qmax=100, Re_ylim=(0, 4e3), tslice=365+174)

## Case 01: Flat topo, KAN_L forcing
cases = [1, 2, 3, 4, 5]
pattern = '../glads/01_kan_forcing/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = 'figures/aux/01_Re_seasonal_1Sept.png'
fig_01 = plot_Re(fnames, figname, Qmin=1, Qmax=100, Re_ylim=(0, 4e3), tslice=365+244)

## Case 01a: KAN increased forcing
cases = [1, 2, 3, 4, 5]
pattern = '../glads/01a_kan_adj_forcing/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = 'figures/aux/01a_Re_seasonal.png'
fig_01 = plot_Re(fnames, figname, Qmin=10, Qmax=200, Re_ylim=(0, 8e3))

## Case 01b: Flat topo, KAN_L forcing, reduced ev
cases = [1, 2, 3, 4, 5]
pattern = '../glads/01b_kan_forcing_ev/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname ='figures/aux/01b_Re_seasonal_ev.png'
fig_01 = plot_Re(fnames, figname, Qmin=10, Qmax=200, Re_ylim=(0, 8e3))

## Case 02a: Synthetic forcing, trough topo
cases = [1, 2, 3, 4, 5]
pattern = '../glads/02a_synth_forcing_trough/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = 'figures/aux/02a_Re_seasonal_trough.png'
fig_02 = plot_Re(fnames, figname, Qmin=1, Qmax=100, Re_ylim=(0, 4e3))

## Case 02b: Synthetic forcing, valley topo
cases = [1, 2, 3, 4, 5]
pattern = '../glads/02b_synth_forcing_valley/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = 'figures/aux/02b_Re_seasonal_valley.png'
fig_02 = plot_Re(fnames, figname, Qmin=1, Qmax=100, Re_ylim=(0, 4e3))

## Case 03a: KAN forcing, trough topo
cases = [1, 2, 3, 4, 5]
pattern = '../glads/03a_kan_forcing_trough/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = 'figures/aux/03a_Re_seasonal_trough.png'
fig_03 = plot_Re(fnames, figname, Qmin=1, Re_ylim=(0, 4e3))

## Case 03b: KAN forcing, valley topo
cases = [1, 2, 3, 4, 5]
pattern = '../glads/03b_kan_forcing_valley/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = 'figures/aux/03b_Re_seasonal_valley.png'
fig_03 = plot_Re(fnames, figname, Qmin=1, Re_ylim=(0, 4e3))

## Case 03c: KAN forcing, trough topo
cases = [1, 2, 3, 4, 5]
pattern = '../glads/03c_kan_forcing_trough2/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = 'figures/aux/03c_Re_seasonal_trough.png'
fig_03 = plot_Re(fnames, figname, Qmin=1, Re_ylim=(0, 4e3))


## Case S01b: Parameter sensitivity
cases = [1, 2, 3, 4, 5]
pattern = '../glads/S01b_parameter_sensitivity/RUN/output_%03d_seasonal.nc'
fnames = [pattern % caseid for caseid in cases]
figname = 'figures/supplement/S01b_Re_parameter.png'
fig_S01 = plot_Re(fnames, figname, Qmin=1, Re_ylim=(0, 4e3), tslice=365+174)

# plt.show()

