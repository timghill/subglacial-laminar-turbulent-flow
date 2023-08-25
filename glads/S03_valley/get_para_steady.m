function para = get_para_steady(config)
% para = get_para_steady(config)
%
% Set para for steady state run

%% Get defaults and unwrap
% addpath('../')
para = get_para(config);
[pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd, pcm] = unwrap_all_para(para);

%% Time
% pt.end = 20*pp.year;
pt.end   = 100*pp.year;  % end time
pt.out_t = pt.start:5*pp.year:pt.end;

% pt.end = 4*pp.day;
% pt.out_t = pt.start:pp.day:pt.end;

%% Synthetic bed topo
addpath('../data/topo_x_squared_para/')
bed_para = 0.05;
pin.bed_elevation = make_anon_fn('@(xy, time) double(bed_elevation_mountain_glacier(xy, time, bed_para))', bed_para);
pin.ice_thickness = make_anon_fn('@(xy, time) double(ice_thickness_mountain_glacier(xy, time, bed_para))', bed_para);

addpath(genpath('../data/kan_l_melt/'))
pin.source_term_s = make_anon_fn('@(xy, time) double(KAN_dist_steady(time, pin, dmesh));', pin, dmesh);
% pin.source_term_s = make_anon_fn('@(xy, time) double(0*xy(:, 1) + 1.158e-6);');
xy = dmesh.tri.nodes;
pin.source_term_c = make_anon_fn('@(time) double(0*xy(:, 1));', xy);

%% Nondimensionalize and re-wrap
[psp, pst, psmd, psin, mesh] = scale_para(pp, pt, pmd, pin, dmesh, ps);
para = wrap_para(pm, pn, pin, ps, pt, pst, psp, pp, mesh, dmesh, psin, pmd, psmd, pcm);
