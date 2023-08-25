function bed = bed_elevation_para(xy, time)
% out = bed_elevation(xy, time)
%
% return bed elevation. Maintain time dependence since the GlaDS code
% assumes this is a function of time

% Bed elevation parameters
trough_zs = [0, 350];
trough_xs = [0, 100e3];
plateau_elev = 350;
ymid = 12.5e3;
scaling_exp = 0.5;
trough_width = 3e3;

% Compute piecewise linear bed trough elevation
x1 = trough_zs(1) + (trough_zs(2) - trough_zs(1))./(trough_xs(2) - trough_xs(1)).*xy(:, 1);
% x2 = trough_zs(2) + (trough_zs(3) - trough_zs(2))./(trough_xs(3) - trough_xs(2)).*(xy(:, 1) - trough_xs(2));
% x1 = x2;
bed_along_trough = x1;

thick_along_trough = plateau_elev - bed_along_trough;

% Gaussian width scales with trough thickness
w = trough_width*mean(thick_along_trough).^(scaling_exp)/2.355./(thick_along_trough.^scaling_exp);
c = ymid + 0.5*ymid*sin(2*pi*xy(:, 1)/50e3);
bed = plateau_elev - thick_along_trough.*exp(-(xy(:, 2)-c).^2/2./w.^2);

end
