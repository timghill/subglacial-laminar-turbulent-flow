"""

Plot image-derived moulin locations and synthetic design

"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib import tri
from scipy import interpolate
import scipy.stats
import scipy.optimize
from osgeo import ogr
import rasterio
import netCDF4 as nc
import cmocean


#######################################
# Define files and parameters
#######################################

# Filepaths
moulin_file = 'moulins_polarstereo.shp'
moulin_density_file = 'yang2016_density.txt'
arcticdem_file = 'arcticdem_mosaic_500m_v3-0_clipped.tif'
synthetic_moulin_file = '../moulins/moulins_068.txt'
synthetic_catchment_file = '../moulins/catchment_map_068.txt'
synthetic_catchment_centers = '../moulins/catchments_068.txt'
figname = '../../../analysis/figures/supplement/moulin_design.png'

# Plot parameters
figsize=(7, 6)
data_color = 'coral'
synth_color = 'teal'

#######################################
# Read moulin locations
#######################################
moulin_shp = ogr.Open(moulin_file)
layer = moulin_shp.GetLayer(0)
n_moulins = layer.GetFeatureCount()
XY = np.zeros((n_moulins, 2))
for i, feature in enumerate(layer):
    gref = feature.GetGeometryRef()
    x, y = gref.GetPoint_2D()
    XY[i, 0] = x
    XY[i, 1] = y

yang16_z_density = np.loadtxt('yang2016_density.txt')
yang16_z_bins = yang16_z_density[:, 0]
yang16_density = yang16_z_density[:, 1]


#######################################
# Read ArcticDEM data
#######################################
dem_img = rasterio.open(arcticdem_file)
dem = dem_img.read(1)
# print(dem_img.shape)
nr, ncol = dem_img.shape
dem_xx = np.zeros((nr, ncol))
dem_yy = np.zeros((nr, ncol))
for ic in range(ncol):
    for ir in range(nr):
        xy = dem_img.xy(ir, ic)
        dem_xx[ir, ic] = xy[0]
        dem_yy[ir, ic] = xy[1]


# #######################################
# # Interpolate elevation to moulin locations
# #######################################
# x_grid = dem_xx[0, :]
# y_grid = dem_yy[:, 0]
# y_grid = y_grid[::-1]
# dem = dem[::-1, :]
# grid_interpolant = interpolate.RegularGridInterpolator((y_grid, x_grid), dem, method='linear', bounds_error=False, fill_value=None)
# moulin_elevation = grid_interpolant(XY[:, ::-1])

## Compute moulin density
surf_fun = lambda x: 6*(np.sqrt(x  +5e3) - np.sqrt(5e3)) + 390

dmesh = nc.Dataset('../mesh/mesh_04.nc')
nodexy = dmesh['tri/nodes'][:].data.T
nodez = surf_fun(nodexy[:, 0])
connect = dmesh['tri/connect'][:].data.T.astype(int) - 1
moulins = np.loadtxt(synthetic_moulin_file)
moulinx = moulins[:, 1]
mouliny = moulins[:, 2]
moulinz = surf_fun(moulinx)
catchment_labels = np.loadtxt(synthetic_catchment_file)
catchments = np.loadtxt(synthetic_catchment_centers)
z_bins = np.arange(400, 2000, 100)
n_moulin_bins = np.zeros(z_bins.shape)
density_bins = np.zeros(len(z_bins)-1)
tri_surf = surf_fun(dmesh['tri/nodes'][0].data.T)
node_area = dmesh['tri/area_nodes'][:].data.T/1e6
for ii in range(len(z_bins)-1):
    moulin_mask = np.logical_and(moulinz>=z_bins[ii], moulinz<z_bins[ii+1])
    n_moulin_ii = len(moulins[moulin_mask])

    node_mask = np.logical_and(tri_surf>=z_bins[ii], tri_surf<z_bins[ii+1])
    band_area = np.sum(node_area[node_mask])
    density_bins[ii] = n_moulin_ii/band_area

z_bin_center = 0.5*(z_bins[:-1] + z_bins[1:])


#######################################
# Make the figure
#######################################
fig = plt.figure(figsize=figsize)
gs = gridspec.GridSpec(4, 2, left=0.11, right=0.9,
    bottom=0.04, top=0.925, wspace=0.05, hspace=0.05,
    width_ratios=(50, 100), height_ratios=(7, 150, 35, 100))

ax0 = fig.add_subplot(gs[1, 0])
ax1 = fig.add_subplot(gs[1, 1])
ax2 = fig.add_subplot(gs[3, :2])
cax = fig.add_subplot(gs[0, 0])

pc = ax0.contourf(dem_xx/1e3, dem_yy/1e3, dem, cmap='gray', 
    levels=np.arange(0, 2100, 100))
ax0.plot(XY[:, 0]/1e3, XY[:, 1]/1e3, linestyle='', marker='.', markersize=5,
    color=data_color, markeredgecolor='k', markeredgewidth=0.25)
ax0.set_aspect('equal')
ax0.set_xlabel('Easting (km)')
ax0.set_ylabel('Northing (km)')
# ax0.text(0.975, 0.975, 'a', ha='right', va='top', fontweight='bold',
#     transform=ax0.transAxes)
ax0.text(0.025, 0.975, 'a', ha='left', va='top', fontweight='bold',
    transform=ax0.transAxes, color='w')

ax1.plot(yang16_z_bins, yang16_density, color=data_color, label='Yang & Smith (2016)')
ax1.plot(z_bin_center, density_bins, color=synth_color, label='Synthetic')
ax1.yaxis.tick_right()
ax1.yaxis.set_label_position('right')
ax1.set_ylabel('Density (km$^{-2}$)')
ax1.set_xlabel('Elevation (m)')
ax1.grid(linestyle=':', linewidth=0.5)
ax1.set_xlim([350, 2000])
ax1.set_ylim([0, 0.065])
ax1.legend(bbox_to_anchor=(0, 1.02, 1, 0.1), ncol=2, mode='expand',
    frameon=False)
ax1.text(0.01, 0.975, 'b', va='top', fontweight='bold',
    transform=ax1.transAxes)


ax2.set_xlim([0, 100])
ax2.set_ylim([0, 25])
ax2.set_aspect('equal')
ax2.set_xlabel('x (km)')
ax2.set_ylabel('y (km)')

# ax3 = ax2.twiny()
# ax2.set_adjustable('box')
# ax3.set_adjustable('box')

xax2 = ax2.secondary_xaxis(location='top')
surf_inv = lambda z: ((z - 390)/6 + np.sqrt(5e3))**2 - 5e3
print(surf_inv(390))
print(surf_inv(1600))
zticks = [400, 600, 800, 1000, 1200, 1400, 1600, 1800]
x2ticks = [surf_inv(z)/1e3 for z in zticks]
xax2.set_ticks(x2ticks)
xax2.set_xticklabels(zticks)
xax2.set_xlabel('Elevation (m)')
# xax2.set_xlabel("new")

mtri = tri.Triangulation(nodexy[:, 0]/1e3, nodexy[:, 1]/1e3, connect)
# ax2.tripcolor(mtri, nodez, vmin=0, vmax=2000, cmap='gray')
# ax2.tricontour(mtri, nodez, levels=np.arange(0, 2100, 100), color='k')
# ax2.tripcolor(mtri, catchment_labels, cmap=cmocean.cm.ice)
ax2.tripcolor(mtri, catchment_labels, cmap='gray')
ax2.plot(catchments[:, 1]/1e3, catchments[:, 2]/1e3, linestyle='',
    marker='.', markersize=10, color='k', label='Catchments',
    markeredgecolor='w', markeredgewidth=0.25)
ax2.plot(moulinx/1e3, mouliny/1e3, linestyle='',
    marker='.', markersize=10, color=synth_color, label='Moulins',
    markeredgecolor='k', markeredgewidth=0.5)
ax2.set_yticks([0, 12.5, 25])
ax2.legend(loc='lower right', markerscale=1.5)
ax2.text(0.01, 0.975, 'c', va='top', fontweight='bold',
    transform=ax2.transAxes, color='w')

cb = fig.colorbar(pc, cax=cax, orientation='horizontal')
cax.xaxis.tick_top()
cax.xaxis.set_label_position('top')
cb.set_label('Elevation (m)')
cb.set_ticks([0, 500, 1000, 1500, 2000])

# ax3 = ax2.twiny()
# ax3.set_xticks([400, 800, 1200, 1600, 2000])
# ax3.set_xlim([surf_fun(0), surf_fun(100e3)])

# ax2.set_aspect('equal')

dmesh.close()

fig.savefig(figname, dpi=600)

plt.show()