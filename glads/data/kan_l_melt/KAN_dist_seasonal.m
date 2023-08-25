function melt = KAN_dist_seasonal(t, pin, dmesh)
    % KAN_dist_seasonal
    %
    % melt = KAN_moulin_seasonal(time, pin, dmesh, ii_moulin, catchmap, ra)
    %
    % Compute moulin inputs for KAN_L forcing case
    %
    % Returns transient, catchment-dependent moulin inputs for each moulin,
    % using degree-day model and KAN_L AWS temperatures

    T_day = 86400;
    T_year = T_day*365;
    ra = 0;
    ramp = max(0, min(t/T_year/25 - 1, 1));

    %% Compute instantaneous melt rate
    xy = dmesh.tri.nodes;
    z = pin.bed_elevation(xy, 0) + pin.ice_thickness(xy, 0);

    surf_melt = KAN_PDD_melt(t, z);

    diurnal = 1 - ra*sin(2*pi*t/T_day);
    melt = surf_melt.*ramp.*diurnal;
    melt(melt<0) = 0;
end
