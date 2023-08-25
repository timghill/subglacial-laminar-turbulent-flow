"""
Random moulin placement for synthetic subglacial drainage simulations.

Moulins density is scaled according to a melt proxy, H_max - h(x),
so that more moulins are placed at the lowest elevations.

Catchment centers are randomly placed with a space-filling sampling design,
and the moulin locations are defined as the lowest node within each
catchment that are not on the boundary and not within a distance
threshold of existing moulins.

Moulin and catchment information are saved in three text files:

moulins_NNN.txt: Columns (index, x, y) for moulin locations

catchments_NNN.txt: Columns (index, x, y) for catchment centers

catchment_map_NNN.txt: Column (index) giving the index of the nearest
    moulin (between 0 and N-1) for each node in the mesh
    index = -1 indicates the node is not drained by a moulin
"""

import argparse

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import tri
import scipy.optimize
import scipy.integrate
import scipy.stats
import netCDF4 as nc

import pyDOE
import cmocean

def main(n_moulin):
    # some thresholds
    # Restrict moulins so they don't form below x_moulin_min. This represents
    # crevassing and the impact of high surface slopes near the terminus
    x_moulin_min = 5e3
    # Domain boundaries - don't place moulins on the boundary
    x_boundary = np.array([0, 100e3])
    y_boundary = np.array([0, 25e3])
    # Minimum distance theshold
    dist_thresh = 2.5e3

    # Read in triangular mesh
    dmesh = nc.Dataset('../mesh/mesh_04.nc')
    nodes = dmesh['tri/nodes'][:].data.T
    connect = dmesh['tri/connect'][:].data.T - 1

    # Find boundary nodes
    node_bmark = np.zeros(len(nodes))
    node_bmark[nodes[:, 0]<=x_moulin_min] = 1
    node_bmark[nodes[:, 0]==x_boundary[1]] = 1
    node_bmark[nodes[:, 1]==y_boundary[0]] = 1
    node_bmark[nodes[:, 1]==y_boundary[1]] = 1

    # Parameters for surface elevation equation
    a = 100e3
    b = 5e3
    k1 = 6
    z0 = 390

    # Latin hypercube sampling (z-y)
    print('Computing LHS design...')
    X_sample = pyDOE.lhs(2, samples=n_moulin, criterion='cm', iterations=1000)
    print('done')

    fig, ax = plt.subplots()
    ax.scatter(X_sample[:, 0], X_sample[:, 1])
    ax.set_title('LHS sample')

    ## Surface elevations for elevation-dependent density
    Z_node = k1*np.sqrt(nodes[:, 0] + b) - k1*np.sqrt(b) + z0

    ## Empirically derived moulin density
    ppf = lambda x: scipy.stats.norm.ppf(x, loc=1138.25, scale=280.12)

    # Scale elevation samples by ppf function
    YZ_scale = np.zeros(X_sample.shape)
    YZ_scale[:, 0] = ppf(X_sample[:, 0])
    YZ_scale[:, 1] = X_sample[:, 1]*(y_boundary[1] - y_boundary[0])

    # Basic evaluation of elevation density mean, sd from sample
    print('Scaled distribution stats:')
    print(np.mean(YZ_scale[:, 0]), np.std(YZ_scale[:, 0]))

    print('Target distribution stats:')
    print('1138.25', '280.12')

    # Compute x locations for elevations by inverting surface elevation function
    X_scale = np.zeros(YZ_scale.shape)
    X_scale[:, 0] = ((YZ_scale[:, 0] - z0)/k1 + np.sqrt(b))**2 - b
    X_scale[:, 1] = YZ_scale[:, 1]

    # Interpolate onto the triangular mesh
    n_nodes = nodes.shape[0]
    catchment_indices = np.zeros(n_moulin).astype(int)
    for i in range(n_moulin):
        node_dist = np.sqrt( (nodes[:, 0] - X_scale[i, 0])**2 + (nodes[:, 1] - X_scale[i, 1])**2 )
        catchment_indices[i] = np.argmin(node_dist)


    # Compute an (N_node x N_moulin) distance matrix and find closets nodes
    dist_matrix = np.sqrt( (np.vstack(nodes[:, 0]) - nodes[catchment_indices, 0])**2\
            + (np.vstack(nodes[:, 1]) - nodes[catchment_indices, 1])**2)
    ixmin = np.argmin(dist_matrix, axis=1).astype(int)
    ixmin[nodes[:, 0]<5e3] = -1

    # Plot catchments as a tripcolor plot
    f, ax = plt.subplots(figsize=(6, 3))
    mtri = tri.Triangulation(nodes[:, 0]/1e3, nodes[:, 1]/1e3, connect)
    ax.tripcolor(mtri, ixmin, vmin=0, vmax=n_moulin-1, shading='flat', cmap=cmocean.cm.gray)
    ax.set_aspect('equal')
    ax.set_xlim(x_boundary/1e3)
    ax.set_ylim(y_boundary/1e3)

    moulin_nodes = np.zeros(catchment_indices.shape, int)
    for i in range(n_moulin):
        catchment_nodes = np.arange(len(ixmin))[ixmin==i].astype(int)

        x_node = nodes[:, 0][catchment_nodes]
        y_node = nodes[:, 1][catchment_nodes]

        dmat = np.sqrt((x_node - np.vstack(nodes[moulin_nodes, 0]))**2 + (y_node - np.vstack(nodes[moulin_nodes, 1]))**2)
        mindist = np.min(dmat, axis=0)
        minmask = np.zeros(len(catchment_nodes))
        minmask[mindist<dist_thresh] = 1
        near_penalty = 1e30*minmask

        x_node = x_node + 1e30*node_bmark[catchment_nodes] + near_penalty
        # Choose the furthest downstream (= lowest elevation) node within
        # the catchment as the moulin. If there are multiple nodes with the
        # same x position (= same elevation), choose randomly amongst
        # these nodes. Add penalties so moulins are not placed on the
        # boundary and so they are not within the minimum distance
        # threshold of a previously placed moulin
        xmin = x_node.min()
        near_indices = np.where(xmin==x_node)[0]
        if (xmin==x_node).astype(int).sum()>1:
            moulin_ix = np.random.choice(near_indices)
        else:
            moulin_ix = near_indices[0]

        moulin_nodes[i] = catchment_nodes[moulin_ix.astype(int)]

    moulin_mat = np.zeros((n_moulin, 3))
    moulin_mat[:, 0] = moulin_nodes
    moulin_mat[:, 1:] = nodes[moulin_nodes]

    catch_mat = np.zeros((n_moulin, 3))
    catch_mat[:, 0] = catchment_indices
    catch_mat[:, 1:] = nodes[catchment_indices]

    # Save moulin indices
    np.savetxt('moulins_%03d.txt' % n_moulin, moulin_mat, '%d')
    np.savetxt('catchments_%03d.txt' % n_moulin, catch_mat, '%d')
    np.savetxt('catchment_map_%03d.txt' % n_moulin, ixmin, '%d')

    # Plot moulins and catchments
    ax.scatter(nodes[catchment_indices, 0]/1e3, nodes[catchment_indices, 1]/1e3, color=[0, 0.5, 0.5], label='Catchment centers')
    ax.scatter(nodes[moulin_nodes, 0]/1e3, nodes[moulin_nodes, 1]/1e3, color=[0.9961, 0.6445, 0.0], label='Moulins')
    ax.legend(bbox_to_anchor=(0.8, -0.35, 0.2, 0.2))
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_title('Moulins and catchments, n=%d' % n_moulin)
    plt.tight_layout()
    f.savefig('catchments_%03d.png' % n_moulin, dpi=300)
    plt.show()

if __name__=='__main__':
    parser = argparse.ArgumentParser(
        prog='PlaceMoulins',)
    parser.add_argument('n_moulins', type=int)
    args = parser.parse_args()
    main(args.n_moulins)
