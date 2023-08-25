function pa = run_job(id)

set_paths;
addpath(genpath('../data/functions/'))

fname_steady = sprintf('output_%03d_steady.nc', id);
fname_seasonal = sprintf('output_%03d_seasonal.nc', id);

% Fixed parameters
config.k_s = 0.1;	% Turbulent conductivity
config.l_c = 10;
config.n_moulin = 68;
config.creep_const_soft = 0;
config.e_v = 1e-4;

config.k_c = 0.2;
bump_scale = 10;
config.h_bed = bump_scale/20;
config.l_bed = bump_scale;
config.alpha = 3/2;
config.beta = 3/2;
config.omega = 1/2000;

config.mesh_nr = id;

if config.alpha<3 && config.omega==0
    % Compute potential gradient for turbulent conductivity scaling
    p_w_min = 40*910*9.81;
    p_w_max = 1520*910*9.81;
    gradphi = (p_w_max - p_w_min)/100e3;
    omega = 1/2000;
    nu = 1.79e-6;
    h3 = nu/(omega)/config.k_s/gradphi;
    k_s = config.k_s * h3^(1 - config.alpha/3) * gradphi^(2 - 3/2);
    config.k_s = k_s;
end

config

% Case-specific filenames
config.fname_steady = fname_steady;
config.fname_seasonal = fname_seasonal;

% Call GlaDS for each parameter set
para_steady = get_para_steady(config);
para_steady.physical;
pa = para_steady;
run_model(para_steady);


% para_seasonal = get_para_seasonal(config);
% run_model(para_seasonal);
