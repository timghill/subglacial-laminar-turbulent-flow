meshfile = load('mesh_refinement.mat');

for ii=1:10
    ncout = sprintf('mesh_refinement_%02d.nc', ii);
    mesh = meshfile.meshes{ii};
    convert_mat_nc(mesh, ncout)
end
