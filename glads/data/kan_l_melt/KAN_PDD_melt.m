function melt = KAN_PDD_melt(t, z, varargin)
    % KAN_PDD_melt Degree-day melt model for KAN_L AWS temperature
    %
    % melt = KAN_PDD_melt(t, z) computes melt rate for default parameters
    %
    % melt = shmip_PDD_melt(t, z, 'DDF', ddf, 'lr', lr, 'DT', dt, 'ra', ra)
    %    specifies the degree-day factor, lapse rate, temperature offset, and
    %    relative diurnal amplitude. Note that the data was generated with a
    %    lapse rate of -0.005 C/km, so be careful using another value
    %
    % Compute instantaneous distributed melt rate with SHMIP PDD
    % melt parameterization
    %
    % Inputs:
    % t: time (seconds)
    % z: elevation (m asl)
    %
    % Options:
    % DDF: Degree-day factor (m w.e./C/day) [Defulat: 0.01 m w.e./C/day)
    % lr: Lapse rate (C/m) [Default: -0.005 C/m]
    % DT: Temperature offset (C) [Defualt: 0 C]
    % ra: Relative diurnal amplitude (-) [Default: 0]


    %% Constants
    day = 86400;            % One day in seconds
    year = 365*day;        % One year in seconds

    %% Default parameters
    DDF = 0.01/day;       % Degree day factor (m/K/s)
    lr = -0.005;           % Lapse rate (K/m)
    DT = 0;                 % Temperature factor
    ra = 0;			% No diurnal melt variation

    if nargin>0
        for ii=1:2:(nargin-2)
            var = varargin{ii};
            val = varargin{ii+1};
            switch var
            case 'DDF'
                DDF = val/day;
            case 'lr'
                lr = val;
            case 'DT'
                DT = val;
            case 'ra'
                ra = val;
            otherwise
                error('Input %s not recognized', var)
            end
        end
    end

    %% Read AWS data
    KAN_data = importdata('KAN_L_2014_temp_clipped.txt');
    tt_days = KAN_data(:, 1);
    sl_temp = KAN_data(:, 2) + DT;

    % Interpolate to find instantaneous melt rate
    tt_day_resid = mod(t, year)/day;
    temp_interp = interp1(tt_days, sl_temp, tt_day_resid, 'linear', 0);

    % Compute elevations and distributed temperature
    diurnal = 1 - ra*sin(2*pi*t/day);
    temp_dist = z*lr + temp_interp;
    melt = DDF*max(0, temp_dist).*diurnal;
    melt(melt<0) = 0;
end
