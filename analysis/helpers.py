"""

Helpful functions and utilities for GlaDS analysis

"""

import numpy as np

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
    z_avg = np.zeros(xmid.shape)
    for i in range(len(xmid)):
        z_avg[i] = metric(z[tf_mask[:, i]])
    
    return xmid, z_avg


def weighted_width_average(xy, z, weights, dx=1e3, xmin=0, xmax=100e3):
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
    z_avg = np.zeros(xmid.shape)
    for i in range(len(xmid)):
        z_avg[i] = np.sum(z[tf_mask[:, i]]*weights[tf_mask[:, i]])/np.sum(weights[tf_mask[:, i]])
    
    return xmid, z_avg