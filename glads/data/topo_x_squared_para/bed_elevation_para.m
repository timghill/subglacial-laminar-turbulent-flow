function out = bed_elevation_para(xy, time)
% out = bed_elevation(xy, time)
%
% return bed elevation. Maintain time dependence since the GlaDS code
% assumes this is a function of time

% out = 0.005*xy(:, 1);
out = 0*xy(:, 1);
