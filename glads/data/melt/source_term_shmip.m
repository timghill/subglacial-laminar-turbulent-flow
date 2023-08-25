function melt = source_term_shmip(xy, time, pin)
% melt = source_term_shmip(xy, time, pin)
%
% Distributed SHMIP degree-day-model melt
% for cases with no moulins

bed = pin.bed_elevation(xy, time);
thick = pin.ice_thickness(xy, time);
z = bed + thick;

melt = shmip_melt(time);

tstart = 25*365*86400;
tend = 50*365*86400;

ramp = (time - tstart)/(tend-tstart);
ramp = max(0, min(1, ramp));

melt = melt * ramp  + 0.0025/86400/365;

% melt = 0.01/86400 * ones(size(xy(:, 1)));
end
