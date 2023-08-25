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

%% Make domain boundary (not so trivial as box case)
boundary_x_forward = linspace(0, 6e3, 501);
boundary_x_reverse = fliplr(boundary_x_forward);
boundary_x_reverse = boundary_x_reverse(2:end-1);
boundary_x = [boundary_x_forward, boundary_x_reverse];

% Topo helper functions
dz = 900;
min_thick = 40;

surface = @(x, y) 100*(x+200).^(1/4) + 1/60*x - 2e10^(1/4);

pa_bench = 0.05;
eps = 1e-16;
f = @(x, pa) (surface(6e3, 0) - pa*6e3)./(6e3).^2 .* x.^2 + pa*x;
g = @(y) 0.5e-6*abs(y).^3;
h = @(x, pa) (-4.5*x/6e3 + 5);%.*(surface(x,0) - f(x, pa))./(surface(x, 0) - f(x, pa_bench) + eps);
bed_fun = @(x, y, pa) f(x, pa) + g(y).*h(x, pa);

ginv = @(x) (x/0.5e-6).^(1/3);
outline_fun = @(x) real(double(ginv((surface(x, 0) - f(x, pa_bench))./(h(x, pa_bench)))));

upper_outline = outline_fun(boundary_x_forward);
lower_outline = -outline_fun(boundary_x_reverse);
outline = [upper_outline, lower_outline];
boundary_xy = [boundary_x; outline]';


figure
subplot(2, 1, 1)
hold on
surf = surface(boundary_x_forward, 0) + dz + min_thick;
bed = bed_fun(boundary_x_forward, 0, pa_bench) + dz;
title('Elevation and Domain boundary')

grid on

plot(boundary_x_forward, surf, 'Color', 'b');
plot(boundary_x_forward, bed, 'Color', 'k');

subplot(2, 1, 2)
plot(boundary_x, outline)

grid on


print('mountain_glacier_elevation_boundary', '-dpng', '-r600')

bmark = 2*ones(length(outline), 1);
bmark(1) = 1;
bmark(end) = 1;
bmark_edge = 2*ones(length(outline), 1);
bmark_edge(1) = 1;
bmark_edge(end) = 1;

maxareas = [1e3, 1.5e3, 2e3, 2.5e3, 5e3];

% cell array holding all the meshes
meshes = {};
for ii=1:length(maxareas)
    meshes{ii} = make_mesh(boundary_xy, bmark, bmark_edge, maxareas(ii));
    figure;
    ax = gca;
    mii = meshes{ii};
    meshtri = mii.tri;
    mesh_plot_tri(ax, meshtri, 1)
    xlim([0, 6e3])
    ylim([-1e3, 1e3])
    title(sprintf('Mesh with %d nodes', meshtri.n_nodes))
    print(sprintf('mountain_glacier_mesh_%03d', ii), '-dpng', '-r600')
end

save([mfiledir, '/valley_mesh.mat'], 'meshes');

