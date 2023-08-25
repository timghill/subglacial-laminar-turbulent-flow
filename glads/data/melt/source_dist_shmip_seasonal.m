function melt = source_dist_shmip_seasonal(xy, time, pin, dmesh)
% melt = source_dist_shmip_seasonal(xy, time, pin)
%
% Compute moulin inputs for SHMIP synthetic case
%
% Returns constant, catchment-dependent moulin inputs for each moulin,
% using SHMIP melt lapse rate

%% Parameters
day = 86400;            % One day in seconds
year = 31536000;        % One year in seconds
DDF = 0.01/86400;       % Degree day factor (m/K/s)
lr = -0.0075;           % Lapse rate (K/m)
DT = 0;                 % Simulate case "D3"
ra = 0;			% No diurnal melt variation


T_day = 86400;
T_year = 86400*365;
ramp = max(0, min(time/T_year/25 - 1, 1));

%% Compute instantaneous melt rate
xy = dmesh.tri.nodes;
z = pin.bed_elevation(xy, 0) + pin.ice_thickness(xy, 0);

temp = -16*cos(2*pi*time/year) - 5 + DT;
surf_melt = DDF*(z*lr + temp);
surf_melt(surf_melt<0) = 0;

melt = surf_melt + 0.01/86400/365;

