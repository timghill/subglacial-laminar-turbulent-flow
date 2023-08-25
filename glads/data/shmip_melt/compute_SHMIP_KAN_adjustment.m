addpath('../kan_l_melt/')

% Set time domain
t0 = 0;
t1 = 365*86400;
dt = 86400;

% Spatial domain
meshes = load('../mesh/mesh.mat');
mesh = meshes.meshes{4};
nodes = mesh.tri.nodes;
areas = mesh.tri.area_nodes;

% Surface elevation
addpath('../topo_x_squared_para/')
z = bed_elevation_flat(nodes,0) + ice_thickness_flat(nodes, 0);

t=0;
KAN_melt = 0;
tt = t0:dt:t1;
while t<t1
    kmelt = KAN_PDD_melt(t, z);
    KAN_melt = KAN_melt + sum(kmelt.*areas)*dt;
    t = t + dt;
end

SHMIP_melt = shmip_total_melt(16, dt, z, areas);
init_ratio = SHMIP_melt/KAN_melt

% Really simple way to estimate best a value
func = @(a) shmip_total_melt(a, dt, z, areas) - KAN_melt;
a0 = 16;
a_estimate = fzero(func, a0)

function melt = shmip_total_melt(a, dt, z, areas)
    T_year = 86400*365;
    ddf = 0.01/86400;
    shmip_lr = 0.0075;
    kan_lr = 0.005;
    dz = 390;

    DT = dz*shmip_lr;
    const = a*(DT - 5)/16;
    T_sl = @(t) -a*cos(2*pi*t/T_year) + const;
    max_temp = T_sl(182.5*86400);
    t=0;
    melt = 0;
    while t<T_year
        smelt = max(0, ddf*(T_sl(t) - z*kan_lr));
	melt = melt + sum(smelt.*areas)*dt;
	t = t + dt;
    end
end
