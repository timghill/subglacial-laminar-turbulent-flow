function melt = source_moulin_shmip_steady(time, pin, dmesh, ii_moulin, catchmap)
    % melt = source_moulin_shmip_steady(time, pin, dmesh, ii_moulin, catchmap)
    % compute steady moulin inputs using SHMIP melt parameterization
    %
    % Uses 25 year winter steady-state + 25 year linear ramp-up of
    % melt intensity to ensure model stability.
    %
    % Uses average melt rate from compute_average_melt
    %
    % See also compute_average_melt

    moulins = zeros(dmesh.tri.n_nodes, 1);
    moulins(ii_moulin) = 1;

    ramp = max(0, min(time/86400/365/25 - 1, 1));

    xy = dmesh.tri.nodes;
    z = pin.bed_elevation(xy, 0) + pin.ice_thickness(xy, 0);

    % Read steady surface melt
    steady_melt = readmatrix('SHMIP_mean_melt.txt');

    area = dmesh.tri.area_nodes;
    catch_melt = integrate_melt_by_catchment(ii_moulin, catchmap, area, steady_melt);

    melt = catch_melt.*ramp;
end
