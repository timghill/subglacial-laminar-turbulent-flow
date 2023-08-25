function bed = bed_elevation_valley(xy, time)
% out = bed_elevation(xy, time)
%
% return bed elevation. Maintain time dependence since the GlaDS code
% assumes this is a function of time

% Bed elevation parameters
trough_dz = -200;
ridge_dz = 400;
yc = 12.5e3;
const_bed = 350;

x = xy(:, 1);
y = xy(:, 2);

bed_ramp = trough_dz*((x - min(x) + 5e3).^0.5 - (5e3).^0.5)./(max(x) - min(x)).^0.5;

ridge_ramp = ridge_dz*((x - min(x) + 5e3).^(1/3) - (5e3).^(1/3))./(max(x) - min(x)).^(1/3);
bed_trough = ridge_ramp.*(y - yc).^2./(max(abs(y - yc))).^2;

bed = const_bed + bed_ramp + bed_trough;

end
