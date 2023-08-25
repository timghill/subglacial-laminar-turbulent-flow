addpath('../shmip_melt')

KAN_data = importdata('KAN_L_2014_temp_clipped.txt');
tt_days = KAN_data(:, 1);
sl_temp = KAN_data(:, 2);

t1 = 0;
t2 = 365*86400;
dt = 3600;
tt = t1:dt:t2;

% Spatial domain
meshes = load('../mesh/mesh.mat');
mesh = meshes.meshes{4};
nodes = mesh.tri.nodes;
areas = mesh.tri.area_nodes;

% Surface elevation
addpath('../topo_x_squared_para/')
z = bed_elevation_flat(nodes,0) + ice_thickness_flat(nodes, 0);

SHMIP_melt = 0;
for ii=1:length(tt)
    t = tt(ii);
    shmip_melt = shmip_PDD_melt(t, z);
    SHMIP_melt = SHMIP_melt + sum(areas.*shmip_melt)*dt;
end

SHMIP_melt

kan_adj_total_melt(0, 86400, z, areas)
% First guess
obj_fun = @(x) kan_adj_total_melt(x, 86400, z, areas) - SHMIP_melt;
melt_factor_guess = 0.5;
melt_factor = fzero(obj_fun, melt_factor_guess);

% Write the adjusted temperature series
adjusted_sl_temp = sl_temp + melt_factor;
mat_write = [tt_days, adjusted_sl_temp]

writematrix(mat_write, 'KAN_L_2014_temp_adjusted.txt')

function kan_melt = kan_adj_total_melt(melt_factor, dt, z, areas)
    t0 = 0;
    t1 = 365*86400;
    tt = t0:dt:t1;

    kan_melt = 0;
    for ii=1:length(tt)
        t = tt(ii);
        kan_melt_ii = KAN_PDD_melt(t, z, 'lr', -0.0075, 'DT', melt_factor);
        kan_melt = kan_melt + sum(areas.*kan_melt_ii)*dt;
    end
end
