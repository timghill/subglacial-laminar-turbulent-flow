%% Compute time-averaged melt using SHMIP PDD model with KAN AWS temperature

ra = 0;         % Diurnal amplitude (0 for no diurnal amplitude)
lr = -0.005;
pa_bench = 0.05;
mesh_nr = 4;    % Mesh index
output_fname = 'KAN_mountain_glacier_mean_melt.txt';  

% Get topo functions
topo_path = '../topo_x_squared_para/';
addpath(genpath(topo_path))

% Load mesh
mesh_struct = load('../mesh/valley_mesh.mat');
meshes = mesh_struct.meshes;
dmesh = meshes{mesh_nr};
xy = dmesh.tri.nodes;

z = bed_elevation_mountain_glacier(xy, 0, pa_bench) + ice_thickness_mountain_glacier(xy, 0, pa_bench);

figure
scatter(xy(:, 1), z)
day = 86400;
year = 365*day;

% Set duration and increment for time integration
t0 = 0;
t1 = year;
dt = day;

% Integration: rectangular midpoint rule
t = t0;
t_melt = 0;
net_melt = 0*xy(:, 1);      % Keep track of total melt
t_melt = 0;                 % Keep track of duration of melt
while t<t1
    inst_melt = KAN_PDD_melt(t, z, 'ra', ra);

    % If any melt happens, increment melt season length
    if max(inst_melt)>0
        t_melt = t_melt + dt;
    end

    net_melt = net_melt + dt*inst_melt;
    t = t + dt;
end

t_melt/86400
mean_melt = net_melt./t_melt;

writematrix(mean_melt, output_fname, 'Delimiter', ',')

% Plot mean melt
[x_sorted, I_sort] = sort(xy(:, 1));
melt_sorted = mean_melt(I_sort);

figure
scatter(xy(:, 1)/1e3, mean_melt*86400)
grid on
xlabel('x (km)')
ylabel('Melt (m w.e. day^{-1})')
title(output_fname, 'Interpreter', 'none')
