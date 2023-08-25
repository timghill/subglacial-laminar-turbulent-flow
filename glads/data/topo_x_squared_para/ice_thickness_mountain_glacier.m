function thick = ice_thickness_mountain_glacier(xy, time, pa)
% surface = ice_thickness_sqrt_shape(xy, time)
%
% Returns the ice thickness. Maintain time dependence since the GlaDS code
% assumes this is a function of time

dz = 900;
min_thick = 40;

surface = @(x, y) 100*(x+200).^(1/4) + 1/60*x - 2e10^(1/4);

surf = surface(xy(:, 1), xy(:, 2)) + dz + min_thick;
bed = bed_elevation_mountain_glacier(xy, time, pa);
thick = surf - bed;

end
