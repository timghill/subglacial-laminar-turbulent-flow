function melt=KAN_moulin_steady(time, pin, dmesh, ii_moulin, catchmap, ra)
    % melt = KAN_moulin_steady(time, pin, dmesh, ii_moulin, catchmap, ra)
    % computes steady moulin inputs based on KAN_L AWS temperatures.
    %
    % Returns constant moulin inputs for each specified moulin,
    % using average KAN_L melt
    %
    % Uses average melt rate from KAN_compute_average_melt
    %
    % See also KAN_compute_average_melt

    ramp = max(0, min(time/86400/365/25 - 1, 1));

    steady_melt = readmatrix('KAN_mean_melt.txt');

    area = dmesh.tri.area_nodes;
    catch_melt = integrate_melt_by_catchment(ii_moulin, catchmap, area, steady_melt);

    melt = catch_melt.*ramp;
end

