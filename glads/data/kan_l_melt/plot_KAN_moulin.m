% function melt=KAN_moulin(time, pin, dmesh, ii_moulin, catchmap, ra)

moulins = load('../moulins/moulins_068.txt');
ii_moulin = moulins(:, 1);
catchmap = load('../moulins/catchment_map_068.txt');
meshes = load('../mesh/mesh.mat');
dmesh = meshes.meshes{4};

addpath('../topo_x_squared_para/')
addpath('../functions/')
%pin.bed_elevation = @(xy, tt) 0*xy(:, 1) + 350;
%pin.ice_thickness = @(xy, tt) 6*(sqrt(xy(:, 1) + 5000) - sqrt(5000)) + 390 - 350;

pin.bed_elevation = @(xy, t) bed_elevation_flat(xy, t);
pin.ice_thickness = @(xy, t) ice_thickness_flat(xy, t);

tt = 86400*(0:0.25:365) + 100*365*86400;
melt = zeros(dmesh.tri.n_nodes, length(tt));
for ii=1:length(tt)
    melt(:, ii) = KAN_moulin_seasonal(tt(ii), pin, dmesh, ii_moulin, catchmap);
end

dt= tt(2) - tt(1);
total_minputs = sum(melt(:))
area = 100e3 * 25e3;
mean_melt = total_minputs/area * dt


t_plt = tt/86400 - tt(1)/86400;
m_plt = 1:length(ii_moulin);
[xx, yy] = meshgrid(t_plt, m_plt);

figure
pcolor(xx, yy, melt(ii_moulin, :))
% caxis([0, 40])
shading flat
cb = colorbar;
cb.Label.String = 'Moulin input (m^3 s^{-1})';
xlabel('Day of year')
ylabel('Moulin')
%cmocean('amp')
colormap(flipud(bone))
print('KAN_moulin_inputs', '-dpng', '-r600')

max(melt(ii_moulin, :), [], 2)

figure
plot(t_plt, melt(ii_moulin(17), :))
