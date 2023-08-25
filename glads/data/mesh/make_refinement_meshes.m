function meshes = make_meshes()
% Make GlaDS meshes
%
% meshes = make_meshes()
%
% Dirichlet on front boundary, no-flux on other boundaries. 25 x 100 km
% rectangular domain.
%
% Specify areas with 'maxareas' variable inside function.

% the path of this mfile (used for paths below)
mfiledir = [fileparts(mfilename('fullpath')), '/'];

% make the box and write it to a file which the meshing library (Triangle) knows
boundary_xy = [0, 0;
               100e3, 0;
               100e3, 25e3;
               0, 25e3];

% boundary marks>0 on edge:
bmark = [1;2;2;1];         % just a mark which is given to the nodes on the boundary
bmark_edge = [2;2;2;1];  % just a mark which is given to the edges on the boundary

% maxareas = 250000000*[0.1, 0.01, 5e-3, 2.5e-3]; % area of triangle , 1e-4
maxareas = 250000000*logspace(-1, log10(5e-4), 10)


% cell array holding all the meshes
meshes = {};
for ii=1:length(maxareas)
    meshes{ii} = make_mesh(boundary_xy, bmark, bmark_edge, maxareas(ii));
end

save([mfiledir, '/mesh_refinement.mat'], 'meshes');

