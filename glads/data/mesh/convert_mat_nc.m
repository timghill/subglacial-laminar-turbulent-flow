function convert_mat_nc(dmesh, ncout)
% function convert_mat_nc(mesh)
%
% Convert a mesh from a .mat file to a
% .nc file for reading into other languages

ncfmt = 'netCDF4';


% Dimensions
ndim = 2;
nodes = dmesh.tri.n_nodes;
elements = dmesh.tri.n_elements;
edges = dmesh.tri.n_edges;

%% File structure
% Groups
%   tri
%   vor
%   perm
% x_extent
% y_extent

nccreate(ncout, 'x_extent', 'Dimensions', {'ndim', ndim}, 'Format', ncfmt)
nccreate(ncout, 'y_extent', 'Dimensions', {'ndim', ndim})

ncwrite(ncout, 'x_extent', dmesh.x_extent)
ncwrite(ncout, 'y_extent', dmesh.y_extent)

%% Group tri
nccreate(ncout, 'tri/type', 'DataType', 'string')
ncwrite(ncout, 'tri/type', dmesh.tri.type)

nccreate(ncout, 'tri/nodes', 'Dimensions', {'nodes', nodes, 'ndim', ndim})
ncwrite(ncout, 'tri/nodes', dmesh.tri.nodes)

nccreate(ncout, 'tri/bmark', 'Dimensions', {'nodes', nodes})
ncwrite(ncout, 'tri/bmark', dmesh.tri.bmark)

nccreate(ncout, 'tri/connect', 'Dimensions', {'elements', elements, 'tri', 3})
ncwrite(ncout, 'tri/connect', dmesh.tri.connect)

nccreate(ncout, 'tri/connect_edge', 'Dimensions', {'edges', edges, 'ndim', 2})
ncwrite(ncout, 'tri/connect_edge', dmesh.tri.connect_edge)

nccreate(ncout, 'tri/bmark_edge', 'Dimensions', {'edges', edges})
ncwrite(ncout, 'tri/bmark_edge', dmesh.tri.bmark_edge)

nccreate(ncout, 'tri/area', 'Dimensions', {'elements'})
ncwrite(ncout, 'tri/area', dmesh.tri.area)

nccreate(ncout, 'tri/edge_length', 'Dimensions', {'edges', edges})
ncwrite(ncout, 'tri/edge_length', dmesh.tri.edge_length)

nccreate(ncout, 'tri/edge_midpoints', 'Dimensions', {'edges', edges, 'ndim', ndim})
ncwrite(ncout, 'tri/edge_midpoints', dmesh.tri.edge_midpoints)

% nccreate(ncout, 'tri/connect_edge_inv', 'Dimensions')

nccreate(ncout, 'tri/area_nodes', 'Dimensions', {'nodes', nodes})
ncwrite(ncout, 'tri/area_nodes', dmesh.tri.area_nodes)

nccreate(ncout, 'tri/connect_edge_el', 'Dimensions', {'edges',edges, 'ndim', ndim})
ncwrite(ncout, 'tri/connect_edge_el', dmesh.tri.connect_edge_el)

% nccreate(ncout, 'tri/neigh_node', 'Dimensions', {'nodes', 'nodes'})
% ncwrite(ncout, 'tri/neigh_node', dmesh.tri.neigh_node)

% nccreate(ncout, 'tri/neigh_edge_node', 'Dimensions', {'edges', edges, 'nodes', nodes})
% ncwrite(ncout, 'tri/neigh_edge_node', dmesh.tri.neigh_edge_node)

% %% Group vor
% 
% nccreate(ncout, 'vor/type', 'DataType', 'string')
% ncwrite(ncout, 'vor/type', dmesh.vor.type)
% 
% nccreate(ncout, 'vor/nodes', 'Dimensions', {})
