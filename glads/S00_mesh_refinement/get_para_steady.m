function para = get_para_steady(config)
% para = get_para_steady(config)
%
% Set para for steady state run

%% Get defaults and unwrap
para = get_para(config);
[pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd, pcm] = unwrap_all_para(para);

%% Time
% pt.end = 20*pp.year;
pt.end   = 100*pp.year;  % end time
pt.out_t = pt.start:5*pp.year:pt.end;

% pt.end = pp.day;
% pt.out_t = 0:43200:pt.end;

%% Synthetic bed topo
addpath('../data/topo_x_squared_para/')
pin.bed_elevation = make_anon_fn('@(xy, time) double(bed_elevation_flat(xy, time))');
pin.ice_thickness = make_anon_fn('@(xy, time) double(ice_thickness_flat(xy, time))');

%% Source functions
n_moulin = config.n_moulin;

addpath(genpath('./data/shmip_melt/'))
pin.source_term_s = make_anon_fn('@(xy, time) double(0.01/86400/365 + 0*xy(:, 1));');
pin.source_term_c = make_anon_fn('@(time) double(source_mesh_refinement(time, pin, dmesh));', pin, dmesh);

% pin.source_term_s = make_anon_fn('@(xy, time) double(source_mesh_refinement(time, xy, pin));', pin);


%% Nondimensionalize and re-wrap
[psp, pst, psmd, psin, mesh] = scale_para(pp, pt, pmd, pin, dmesh, ps);
para = wrap_para(pm, pn, pin, ps, pt, pst, psp, pp, mesh, dmesh, psin, pmd, psmd, pcm);
