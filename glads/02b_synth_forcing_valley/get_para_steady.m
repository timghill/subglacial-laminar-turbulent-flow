function para = get_para_steady(config)
% para = get_para_steady(config)
%
% Set para for steady state run

%% Get defaults and unwrap
addpath('../')
para = get_para(config);
[pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd, pcm] = unwrap_all_para(para);

%% Time
% pt.end = 20*pp.year;
pt.end   = 100*pp.year;  % end time
pt.out_t = pt.start:5*pp.year:pt.end;

%% Synthetic bed topo
addpath('../data/topo_x_squared_para/')
pin.bed_elevation = make_anon_fn('@(xy, time) double(bed_elevation_valley(xy, time))');
pin.ice_thickness = make_anon_fn('@(xy, time) double(ice_thickness_valley(xy, time, pin))', pin);

op = {'freeze-on', 'freeze-melt', 'none'};;
pp.flags.include_sheet_diss_press = op{1}; % if true, allows freezing shut of the sheet on reverse slopes

%% Source functions
n_moulin = config.n_moulin;
moulindata = readmatrix(sprintf('../data/moulins/moulins_%03d.txt', n_moulin));
catchmap = readmatrix(sprintf('../data/moulins/catchment_map_%03d.txt', n_moulin));
ii_moulin = moulindata(:, 1) + 1;

addpath(genpath('../data/shmip_melt/'))
pin.source_term_s = make_anon_fn('@(xy, time) double(0.01/86400/365 + 0*xy(:, 1));');
pin.source_term_c = make_anon_fn('@(time) double(source_moulin_shmip_adj_steady(time, pin, dmesh, ii_moulin, catchmap));', pin, dmesh, ii_moulin, catchmap);

pp.float_frac = 1;
pin.bc_dirichlet_phi = make_anon_fn('@(xy, time, bmark, phi_0, phi_m) double(pp.float_frac * (phi_0-phi_m) + phi_m)', pp);



%% Nondimensionalize and re-wrap
[psp, pst, psmd, psin, mesh] = scale_para(pp, pt, pmd, pin, dmesh, ps);
para = wrap_para(pm, pn, pin, ps, pt, pst, psp, pp, mesh, dmesh, psin, pmd, psmd, pcm);
