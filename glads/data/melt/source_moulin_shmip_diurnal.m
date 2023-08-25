function melt = source_moulin_shmip_diurnal(time, pin, dmesh, ii_moulin, catchmap, ra)
% melt = source_moulin_shmip(xy, time, pin)
%
% Compute moulin inputs for SHMIP synthetic case
%
% Returns constant, catchment-dependent moulin inputs for each moulin,
% using SHMIP melt lapse rate

moulins = zeros(dmesh.tri.n_nodes, 1);
moulins(ii_moulin) = 1;

T_day = 86400;
T_year = 86400*365;
ramp = max(0, min(time/T_year/25 - 1, 1));

xy = dmesh.tri.nodes;
z = pin.bed_elevation(xy, 0) + pin.ice_thickness(xy, 0);

surf_melt = shmip_melt(time);   % Constant, melt season-averaged melt

area = dmesh.tri.area_nodes;
catch_melt = zeros(size(moulins));

for i=1:length(ii_moulin)
    mask = catchmap==(i-1);
    catch_area = area(mask);
    catch_melt(ii_moulin(i)) = sum(catch_area.*surf_melt(mask));
end

diurnal = 1 - ra*sin(2*pi*time/T_day);
melt = catch_melt.*ramp.*diurnal;
melt(melt<0) = 0;
end
