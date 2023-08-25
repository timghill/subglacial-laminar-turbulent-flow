% Compare SHMIP-adjusted and KAN melt forcing

addpath(genpath('../kan_l_melt/'))

addpath('../topo_x_squared_para/')

addpath('../functions/')

pin.bed_elevation = @(xy, t) bed_elevation_flat(xy, t);
pin.ice_thickness = @(xy, t) ice_thickness_flat(xy, t);

meshes = load('../mesh/mesh.mat');
dmesh = meshes.meshes{4};

n_moulin = 68;
moulindata = readmatrix(sprintf('../moulins/moulins_%03d.txt', n_moulin));
catchmap = readmatrix(sprintf('../moulins/catchment_map_%03d.txt', n_moulin));
ii_moulin = moulindata(:, 1) + 1;


dt = 86400;
tt = (0:dt:(365*86400)) + 86400*365*100;
KAN_moulins = 0;
SHMIP_moulins = 0;
for ii=1:length(tt)
    ti = tt(ii);

    KAN_moulins = KAN_moulins + dt*KAN_moulin_seasonal(ti, pin, dmesh, ii_moulin, catchmap);

    SHMIP_moulins = SHMIP_moulins + dt*source_moulin_shmip_adj_seasonal(ti, pin, dmesh, ii_moulin, catchmap);

end

KAN_moulins(ii_moulin)
SHMIP_moulins(ii_moulin)

scatter(KAN_moulins, SHMIP_moulins)
hold on
grid on
plot(KAN_moulins, KAN_moulins, 'k')
xlabel('KAN')
ylabel('SHMIP adjusted')
