function pa = run_job(k_c, bump_scale, alpha, beta, omega, id)

set_paths;
addpath(genpath('../data/functions/'))

fname_steady = sprintf('output_%03d_steady.nc', id);
fname_seasonal = sprintf('output_%03d_seasonal.nc', id);

% Fixed parameters
% config.k_s = 0.1;	% Laminar conductivity
% config.k_s = 0.001;
config.k_s = 1e-4;
config.l_c = 2;
config.n_moulin = 68;
config.creep_const_soft = 0;
config.mesh_nr = 4;

% Tuning parameters
config.k_c = k_c;
config.h_bed = bump_scale/20;
config.l_bed = bump_scale;
config.alpha = alpha;
config.beta = beta;
config.omega = omega;
config.e_v = 1e-4;

if config.alpha<3 && config.omega==0
    % Compute potential gradient for turbulent conductivity scaling
    phi_min = 1000*9.8*900 + 0;
    phi_max = 1000*9.8*1500 + 0;
    gradphi = (phi_max - phi_min)/6e3;
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
run_model(para_steady);

% para_steady
% para_steady.physical.flags

para_seasonal = get_para_seasonal(config);
run_model(para_seasonal);


