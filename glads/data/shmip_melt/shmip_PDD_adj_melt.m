function melt = shmip_PDD_melt(t, z, varargin)
    % SHMIP_PDD_MELT SHMIP degree-day melt model
    % 
    % melt = shmip_PDD_melt(t, z) computes melt rate for default parameters
    %
    % melt = shmip_PDD_melt(t, z, 'DDF', ddf, 'lr', lr, 'DT', dt, 'ra', ra)
    %    specifies the degree-day factor, lapse rate, temperature offset, and
    %    relative diurnal amplitude
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
    % lr: Lapse rate (C/m) [Default: -0.0075 C/m]
    % DT: Temperature offset (C) [Defualt: 0 C]
    % ra: Relative diurnal amplitude (-) [Default: 0]

    %% Constants
    day = 86400;            % One day in seconds
    year = 365*day;        % One year in seconds
    
    %% Default parameters
    DDF = 0.01/day;       % Degree day factor (m/K/s)
    % kan_lr = -0.005;           % Lapse rate (K/m)
    kan_lr = -0.005;
    shmip_lr = -0.0075;
    DT = -shmip_lr*390;                 % Simulate case "D3"
    ra = 0;			% No diurnal melt variation
    a = 9.0684;
    % a = 16;

    if nargin>0
        for ii=1:2:(nargin-2)
            var = varargin{ii};
            val = varargin{ii+1};
            switch var
            case 'DDF'
                DDF = val/day;
            case 'lr'
                lr = val;
         %   case 'DT'
         %       DT = val;
            case 'ra'
                ra = val;
            otherwise
                error('Input %s not recognized', var)
            end
        end
    end

    %% Compute instantaneous melt rate

    % Temperature parameterization
    % temp = -16*cos(2*pi*t/year) - 5 + DT;
    temp = -a*cos(2*pi*t/year) + a*(DT - 5)/16;

    % Degree-day model with lapse rate and diurnal cycle
    melt = DDF*(z*kan_lr + temp);
    diurnal = 1 - ra*sin(2*pi*t/day);

    melt = melt.*diurnal;
    melt(melt<0) = 0;
    end
