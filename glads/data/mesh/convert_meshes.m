meshfile = load('mesh.mat');

for ii=1:4
    ncout = sprintf('mesh_%02d.nc', ii);
    mesh = meshfile.meshes{ii};
    convert_mat_nc(mesh, ncout)
end