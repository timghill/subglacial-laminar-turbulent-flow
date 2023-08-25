function melt = source_moulin_shmip(time, pin, dmesh, ii_moulin, catchmap)
% melt = source_moulin_shmip(xy, time, pin)
%
% Compute moulin inputs for SHMIP synthetic case
%
% FIRST
% Returns constant and prescribed moulin inputs for each moulin with
% a 1 year ramp for stability
%
% SECOND
% Returns constant, catchment-dependent moulin inputs for each moulin,
% using SHMIP melt lapse rate

%% EASY WAY
% ramp = max(0, min(time/86400/365, 1));
% m_const = 25;
% melt = ramp*m_const*moulins;

%% PROPER WAY
moulins = zeros(dmesh.tri.n_nodes, 1);
moulins(ii_moulin) = 1;

ramp = max(0, min(time/86400/365/25 - 1, 1));

xy = dmesh.tri.nodes;
z = pin.bed_elevation(xy, 0) + pin.ice_thickness(xy, 0);

surf_melt = shmip_melt_annual(time);

area = dmesh.tri.area_nodes;
catch_melt = zeros(size(moulins));

for i=1:length(ii_moulin)
    mask = catchmap==(i-1);
    catch_area = area(mask);
    catch_melt(ii_moulin(i)) = sum(catch_area.*surf_melt(mask));
end

melt = catch_melt.*ramp;
end
