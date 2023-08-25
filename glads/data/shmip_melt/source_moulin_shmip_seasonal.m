function melt = source_moulin_shmip_seasonal(t, pin, dmesh, ii_moulin, catchmap, ra)
    % SOURCE_MOULIN_SHMIP_SEASONAL
    %
    % melt = source_moulin_shmip_seasonal(time, pin, dmesh, ii_moulin, catchmap, ra)
    %
    % Compute moulin inputs for SHMIP synthetic case
    %
    % Returns transient, catchment-dependent moulin inputs for each moulin,
    % using SHMIP melt lapse rate

    %% Parameters
    ra = 0;                 % No diurnal melt variation
    lr = -0.0075;           % Lapse rate (C/m)

    T_day = 86400;
    T_year = T_day*365;
    ramp = max(0, min(t/T_year/25 - 1, 1));

    %% Compute instantaneous melt rate
    xy = dmesh.tri.nodes;
    z = pin.bed_elevation(xy, 0) + pin.ice_thickness(xy, 0);

    DT = -min(z).*lr;

    surf_melt = shmip_PDD_melt(t, z, 'DT', DT, 'ra', ra);

    %% Compute catchment melt and moulin inputs
    area = dmesh.tri.area_nodes;
    catch_melt = integrate_melt_by_catchment(ii_moulin, catchmap, area, surf_melt);

    diurnal = 1 - ra*sin(2*pi*t/T_day);
    melt = catch_melt.*ramp.*diurnal;
    melt(melt<0) = 0;
end
