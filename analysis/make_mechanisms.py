"""

Call plot_mechanisms.py

"""

from matplotlib import pyplot as plt

from plot_mechanisms import plot_mechanisms
import defaults

## Case 00a: Flat topo, SHMIP forcing
cases = [1, 2, 3, 4, 5]
fnames = ['../glads/00_synth_forcing/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/main/00_mechanisms.png'
models = defaults.labels
fig_00 = plot_mechanisms(fnames, figname, models, omega=1/2000,
    k_turb=0.00954298166213413, tslice=530)

## Case 00a: Flat topo, SHMIP forcing
cases = [1, 2, 3, 4, 5]
fnames = ['../glads/00a_shmip_forcing/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/00a_mechanisms_shmip_forcing.png'
models = defaults.labels
fig_00 = plot_mechanisms(fnames, figname, models, omega=1/2000,
    k_turb=0.00954298166213413, tslice=530)

## Case 00c: marine-terminating
cases = [1, 2, 3, 4, 5]
fnames = ['../glads/00c_synth_marine/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/00c_mechanisms_marine.png'
models = defaults.labels
fig_00 = plot_mechanisms(fnames, figname, models, omega=1/2000,
    k_turb=0.00954298166213413, tslice=530)

## Case 1: KAN forcing
cases = [1, 2, 3, 4, 5]
fnames = ['../glads/01_kan_forcing/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/01_mechanisms.png'
models = defaults.labels
fig_01 = plot_mechanisms(fnames, figname, models, omega=1/2000,
    k_turb=0.00954298166213413, tslice=365+190, labels=models)

cases = [3, 6, 2, 4, 5]
fnames = ['../glads/01_kan_forcing/RUN/output_%03d_seasonal.nc'%caseid for caseid in cases]
figname = 'figures/aux/01_mechanisms_alpha.png'
models = ['Laminar', 'Laminar 4', 'Turbulent 3/2', 'Transition 5/4' ,'Transition 3/2']
fig_01 = plot_mechanisms(fnames, figname, models, omega=1/2000,
    k_turb=0.00954298166213413, tslice=365+190, labels=models)


plt.show()

