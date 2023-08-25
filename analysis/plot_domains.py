"""
Plot bed elevation and ice thickness for bed topo realizations
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import tri
from matplotlib.transforms import Bbox
from matplotlib.gridspec import GridSpec
import netCDF4 as nc
import cmocean
from palettes.code import palettes

# Get mesh info
dmesh = nc.Dataset('../glads/data/mesh/mesh_04.nc')
nodes = dmesh['tri/nodes'][:].data.T
connect = dmesh['tri/connect'][:].data.T.astype(int) - 1

# Get bed elevation from model outputs
out_trough1 = nc.Dataset('../glads/02a_synth_forcing_trough/RUN/output_003_seasonal.nc')
out_trough2 = nc.Dataset('../glads/03c_kan_forcing_trough2/RUN/output_003_seasonal.nc')
out_valley = nc.Dataset('../glads/02b_synth_forcing_valley/RUN/output_003_seasonal.nc')

bed_flat = 350 + 0*nodes[:, 0]
bed_trough1 = out_trough1['bed'][:].data
bed_trough2 = out_trough2['bed'][:].data
bed_valley = out_valley['bed'][:].data

# All cases share surface elevation
surf = 6*(np.sqrt(nodes[:, 0] + 5e3) - np.sqrt(5e3)) + 390

# Plotting
mesh_tri = tri.Triangulation(nodes[:, 0]/1e3, nodes[:, 1]/1e3, connect)
bed_levels = np.arange(0, 550, 50)
thick_levels = np.arange(0, 2100, 100)

# Plot bed elevation and ice thickness
fig2 = plt.figure(figsize=(7, 5))

gs = GridSpec(5, 2, height_ratios=(10, 100, 100, 100, 100),
    left=0.1, bottom=0.1, right=0.97, top=0.9,
    hspace=0.05, wspace=0.075)

axs = np.array([[fig2.add_subplot(gs[i+1,j]) for j in range(2)] for i in range(4)]).T

bed_cmap = palettes.get_cmap('BrownEarth').reversed()
thick_cmap = palettes.get_cmap('BlueIce')

bed_pcolor = axs[0, 0].tricontourf(mesh_tri, bed_flat, cmap=bed_cmap, levels=bed_levels)
thick_pcolor = axs[1, 0].tricontourf(mesh_tri, surf - bed_flat, cmap=thick_cmap,  levels=thick_levels)

axs[0, 1].tricontourf(mesh_tri, bed_trough1, cmap=bed_cmap, levels=bed_levels)
axs[1, 1].tricontourf(mesh_tri, surf - bed_trough1, cmap=thick_cmap, levels=thick_levels)

axs[0, 2].tricontourf(mesh_tri, bed_trough2, cmap=bed_cmap, levels=bed_levels)
axs[1, 2].tricontourf(mesh_tri, surf - bed_trough2, cmap=thick_cmap, levels=thick_levels)

axs[0, 3].tricontourf(mesh_tri, bed_valley, cmap=bed_cmap,  levels=bed_levels)
axs[1, 3].tricontourf(mesh_tri, surf - bed_valley, cmap=thick_cmap, levels=thick_levels)

for ax in axs.flat:
    ax.set_aspect('equal')
    ax.set_xlim([0, 100])
    ax.set_ylim([0, 25])
    ax.set_yticks([0, 12.5, 25])

axT = axs.T
axT[0, 0].set_xticklabels([])
axT[0, 1].set_xticklabels([])
axT[0, 0].text(-0.2, 0.95, 'a', transform=axT[0, 0].transAxes, fontweight='bold')


axT[1, 0].set_xticklabels([])
axT[1, 1].set_xticklabels([])
axT[2, 0].set_xticklabels([])
axT[2, 1].set_xticklabels([])
axT[1, 0].text(-0.2, 0.95, 'b', transform=axT[1, 0].transAxes, fontweight='bold')
# axT[0, 2].set_xticklabels([])

axT[0, 1].set_yticklabels([])
axT[1, 1].set_yticklabels([])
axT[2, 1].set_yticklabels([])
axT[3, 1].set_yticklabels([])
axT[2, 0].text(-0.2, 0.95, 'c', transform=axT[2, 0].transAxes, fontweight='bold')
axT[3, 0].text(-0.2, 0.95, 'd', transform=axT[3, 0].transAxes, fontweight='bold')

axT[0, 0].set_ylabel('y (km)')
axT[1, 0].set_ylabel('y (km)')
axT[2, 0].set_ylabel('y (km)')
axT[3, 0].set_ylabel('y (km)')

axT[-1, 0].set_xlabel('x (km)')
axT[-1, 1].set_xlabel('x (km)')

cax1 = fig2.add_subplot(gs[0, 0])
cax2 = fig2.add_subplot(gs[0, 1])
cbar1 = fig2.colorbar(bed_pcolor, cax=cax1, orientation='horizontal')
cax1.xaxis.tick_top()
cax1.xaxis.set_label_position('top')

cbar2 = fig2.colorbar(thick_pcolor,cax=cax2, orientation='horizontal')
cax2.xaxis.tick_top()
cax2.xaxis.set_label_position('top')

cbar1.set_label('Bed Elevation (m)')
cbar2.set_label('Thickness (m)')
fig2.savefig('figures/supplement/bed_thickness.png', dpi=600)

plt.show()
