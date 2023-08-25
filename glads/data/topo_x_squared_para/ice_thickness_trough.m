function thick = ice_thickness_para(xy, time, pin)
% surface = ice_thickness_sqrt_shape(xy, time)
%
% Returns the ice thickness. Maintain time dependence since the GlaDS code
% assumes this is a function of time

s = 6*(sqrt(xy(:, 1) + 5e3) - sqrt(5e3)) + 390;
b = pin.bed_elevation(xy, time);
thick = s - b;
end
