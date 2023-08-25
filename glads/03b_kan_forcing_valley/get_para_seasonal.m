function para = get_para_seasonal(config)
% para = get_para_diurnal(mesh_nr, n_moulin, omega, kc, ks, alpha, beta, hr, fname_steady, filename)
%
% Set para for diurnal runs

fname_steady = config.fname_steady;
filename = config.fname_seasonal;

%% Get defaults and unwrap
addpath('../')
para = get_para(config);
[pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd, pcm] = unwrap_all_para(para);

pm.model_run_descript = 'Run seasonal';
pm.save_filename = [pm.dir.model_save, filename];

%% Time
pt.start = 100*pp.year;
pt.end   = pt.start + 2*pp.year;  % end time
pt.out_t = pt.start : 1*pp.day : pt.end;

%% Synthetic bed topo
addpath('../data/topo_x_squared_para/')
pin.bed_elevation = make_anon_fn('@(xy, time) double(bed_elevation_valley(xy, time))');
pin.ice_thickness = make_anon_fn('@(xy, time) double(ice_thickness_valley(xy, time, pin))', pin);


%% Source functions
addpath('../data/kan_l_melt/')
n_moulin = config.n_moulin;
moulindata = readmatrix(sprintf('../data/moulins/moulins_%03d.txt', n_moulin));
catchmap = readmatrix(sprintf('../data/moulins/catchment_map_%03d.txt', n_moulin));
ii_moulin = moulindata(:, 1) + 1;

% Distributed basal melt
pin.source_term_s = make_anon_fn('@(xy, time) double(0.01/86400/365 + 0*xy(:, 1));');

% Moulin inputs will need to be adjusted for diurnal simulations
pin.source_term_c = make_anon_fn('@(time) double(KAN_moulin_seasonal(time, pin, dmesh, ii_moulin, catchmap));', pin, dmesh, ii_moulin, catchmap);

%% Numerics
st = {'ode15s', 'ode23s', 'ode23t', 'odebim'};  % can also use ode23t, ode23s but ode15s is probbaly best
pn.ts.ode15s.stepper = st{1};

%% Initial conditions from file
% Since we are saving model runs as *.nc instead of *.mat,
% we need to open the *.nc file and put the relevant
% variables into a struct

IC_nc_fname = [pm.dir.model_save, fname_steady];

IC_phi = ncread(IC_nc_fname, 'phi');
IC_h = ncread(IC_nc_fname, 'h_sheet');
IC_S = ncread(IC_nc_fname, 'S_channel');

% Need to non-dimensionalize these fields!
pa_steady = get_para_steady(config);
ps_steady = pa_steady.scale;

IC_phi = IC_phi./ps_steady.phi;
IC_h = IC_h./ps_steady.h;
IC_S = IC_S./ps_steady.S;

IC_fields = struct;
IC_fields.phi = IC_phi(:, end);
IC_fields.h_sheet = IC_h(:, end);
IC_fields.S_channel = IC_S(:, end);

pm.IC_from_file = true;
pm.file.IC_file = IC_fields;

%% Nondimensionalize and re-wrap
[psp, pst, psmd, psin, mesh] = scale_para(pp, pt, pmd, pin, dmesh, ps);
para = wrap_para(pm, pn, pin, ps, pt, pst, psp, pp, mesh, dmesh, psin, pmd, psmd, pcm);
