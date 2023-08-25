"""

Extend moulin placement scheme to meshes used for mesh refinement

"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import tri
import netCDF4 as nc

## Read reference mesh data and moulins
ref_mesh_file = '../mesh/mesh_04.nc'
ref_mesh = nc.Dataset(ref_mesh_file, 'r')
ref_nodes = ref_mesh['tri/nodes'][:].data.T

ref_moulin_file = 'moulins_068.txt'
ref_moulin_data = np.loadtxt(ref_moulin_file)
ref_moulin_xy = ref_moulin_data[:, 1:3]
ref_moulin_index = ref_moulin_data[:, 0].astype(int)

## Read meshes for mesh refinement
n_meshes = 10
n_moulin = ref_moulin_index.shape[0]
for i in range(n_meshes):
# for i in [1]:
    mesh_file = '../mesh/mesh_refinement_%02d.nc' % (i+1)
    mesh = nc.Dataset(mesh_file, 'r')
    nodes = mesh['tri/nodes'][:].data.T

    # Interpolate reference moulin positions onto each new mesh
    moulin_xy = np.zeros(ref_moulin_xy.shape)
    moulin_ix = np.zeros(ref_moulin_index.shape)
    for j in range(n_moulin):
        dist = np.sqrt( (ref_moulin_xy[j, 0] - nodes[:, 0])**2 + (ref_moulin_xy[j, 1] - nodes[:, 1])**2)

        # Enforce no moulins on boundary
        dist[nodes[:, 0]==0] = 1e6
        dist[nodes[:, 0]==100e3] = 1e6
        dist[nodes[:, 1]==0] = 1e6
        dist[nodes[:, 1]==25e3] = 1e6

        ixmin = np.argmin(dist)
        moulin_xy[j, :] = nodes[ixmin, :]
        moulin_ix[j] = ixmin

    moulin_data = np.zeros((moulin_xy.shape[0], 3))
    moulin_data[:, 0] = moulin_ix
    moulin_data[:, 1:] = moulin_xy
    moulin_file = 'mesh_refinement/moulins_%05d.txt' % nodes.shape[0]
    np.savetxt(moulin_file, moulin_data, '%d')

    fig, axs = plt.subplots(nrows=2, sharex=True)
    ax1,ax2 = axs
    ax1.scatter(ref_nodes[:, 0]/1e3, ref_nodes[:, 1]/1e3, color='gray', s=1, zorder=0)
    ax1.scatter(ref_moulin_xy[:, 0]/1e3, ref_moulin_xy[:, 1]/1e3, zorder=2)
    ax1.set_aspect('equal')
    ax1.set_yticks([0, 12.5, 25])
    ax1.set_xlim([0, 100])
    ax1.set_ylim([0, 25])

    ax2.scatter(nodes[:, 0]/1e3, nodes[:, 1]/1e3, color='gray', s=1, zorder=0)
    ax2.scatter(moulin_xy[:, 0]/1e3, moulin_xy[:, 1]/1e3, zorder=2)
    ax2.set_aspect('equal')
    ax2.set_xlim([0, 100])
    ax2.set_ylim([0, 25])
    ax2.set_yticks([0, 12.5, 25])

    figname = 'moulins_%02d' % (i+1)
    fig.savefig(figname, dpi=600)

    mesh.close()

plt.show()