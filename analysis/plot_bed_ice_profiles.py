"""

Plot bed and ice surface elevation profiles

"""

import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.transforms import Bbox

from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Domain topography
bed_fun = lambda x: 350 + 0*x
surf_fun = lambda x: 6*(np.sqrt(x  +5e3) - np.sqrt(5e3)) + 390

# Graphics
# ice_color = 'lightcyan'
ice_color = 'lightblue'
bed_color = 'peru'
sky_color = 'whitesmoke'
line_color = 'k'
# band_color = 'lightslategrey'
band_color = 'steelblue'
figname = 'figures/main/domain_profile.png'
figsize=(5, 4.5)

bands = [15, 30, 70]
band_width = 5
band_background = [0.5, 0.5, 0.5, 0.25]
# band_color = 'dimgray'

# Compute topo on x profile
x = np.linspace(0, 100e3, 101)
y = np.linspace(0, 25e3, 5)
[xx, yy] = np.meshgrid(x, y)
bed = bed_fun(xx)
surf = surf_fun(xx)

# Make the figure
fig = plt.figure(figsize=figsize)
ax = fig.add_subplot(projection='3d', computed_zorder=False)

# Plot bed
ax.plot_surface(xx/1e3, yy/1e3, bed, color=bed_color, edgecolor='none', zorder=0)
xs = [0, 100]
ys = [0, 0]
zzs = 0*xs
[xxs, yys] = np.meshgrid(xs, ys)
print(xxs)
print(yys)
zzs = np.array([[0, 0], [350, 350]])
ax.plot_surface(xxs, yys, zzs, color=bed_color, zorder=1, antialiased=False, alpha=1)

xxs = np.array([[0, 0], [0, 0]])
yys = np.array([[0, 25], [0, 25]])
zzs = np.array([[0, 0], [350, 350]])
ax.plot_surface(xxs, yys, zzs, color=bed_color, zorder=1, antialiased=False, alpha=1)

# Plot surface
[xx1, yy1] = np.meshgrid(x, [0, 0])
z1 = surf_fun(xx1)
z1[0] = 390
print(xx1)
print(yy1)
print(z1)
ax.plot_surface(xx/1e3, yy/1e3, surf, color=ice_color, edgecolor=ice_color, alpha=1, zorder=0, antialiased=False, linewidth=0.25)
ax.plot_surface(xx1/1e3, yy1/1e3, z1, color=ice_color, edgecolor=ice_color, alpha=1, zorder=0, antialiased=False, linewidth=0.25)

# Plot bands
for xb in bands:
    xs = np.linspace(xb-2.5, xb+2.5, 101)
    ys = np.array([0, 25])
    [xxs, yys] = np.meshgrid(xs, ys)

    zs = np.array([[350, surf_fun(xb*1e3)], [350, surf_fun(xb*1e3)]])

    ax.plot([xb, xb], [0, 25], [surf_fun(xb*1e3), surf_fun(xb*1e3)], color=line_color, zorder=5, linewidth=2)
    ax.plot([xb, xb], [0, 0], [390, surf_fun(xb*1e3)], color=line_color, zorder=5, linewidth=2)
    zzs = surf_fun(xxs*1e3)
    ax.plot_surface(xxs, yys, zzs+1, zorder=2, color=band_color, alpha=1, antialiased=False)
    z3 = zzs.copy()
    z3[0] = 350
    ax.plot_surface(xxs, yys*0, z3, zorder=2, color=band_color, alpha=1, antialiased=False)

# Plot moulins
dmesh = nc.Dataset('../glads/data/mesh/mesh_04.nc')
moulins = np.loadtxt('../glads/data/moulins/moulins_068.txt')
print(moulins)
moulinx = moulins[:, 1]
mouliny = moulins[:, 2]
moulinz = surf_fun(moulinx)

ax.plot(moulinx/1e3, mouliny/1e3, moulinz, marker='.', color='k', zorder=6, linestyle='', markersize=5)


ax.set_xticks([0, 20, 40, 60, 80, 100])
ax.set_yticks([0, 12.5, 25])
ax.set_zticks([0, 500, 1000, 1500, 2000])

ax.set_xlabel('Distance from terminus (km)', labelpad=20)
ax.set_ylabel('Distance (km)')
ax.zaxis.set_rotate_label(False)  # disable automatic rotation
ax.set_zlabel('z (m)', rotation=90)

ax.view_init(elev=24, azim=-125) #Works!
ax.set_box_aspect((4, 1, 1))
ax.set_aspect('equalxy')

ax.set_position(Bbox.from_extents(0.1, 0.315, 1., 1.2))
# ax.text(0.05, 0.95, 'a', fontweight='bold')

## Compute moulin density
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

yang_2016_density = np.loadtxt('../glads/data/yang2016_moulins/yang2016_density.txt')

z_bin_center = 0.5*(z_bins[:-1] + z_bins[1:])
ax2 = fig.add_subplot()
ax2.set_position(Bbox.from_bounds(0.25, 0.125, 0.5, 0.25))
ax2.plot(yang_2016_density[:, 0], yang_2016_density[:, 1], label='Target density', color='k')
ax2.plot(z_bin_center, density_bins, label='Synthetic density', color='k', linestyle='dashed')
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.set_xlim([200, 2000])
ax2.set_ylim([0, 0.07])
ax2.text(-0.25, 2.75, 'a', transform=ax2.transAxes, fontweight='bold')
ax2.text(-0.25, 1, 'b', transform=ax2.transAxes, fontweight='bold')

ax2.legend(bbox_to_anchor=(-0.15, 0.9, 1.3, 0.2), ncol=2, mode='expand', frameon=False)

ax2.set_xlabel('Elevation (m)')
ax2.set_ylabel(r'Density (km$^{-2}$)')
ax2.grid()

fig.savefig(figname, dpi=600)

plt.show()

