function melt = integrate_melt_by_catchment(moulin_indices, catchments, node_areas, catchment_melt)
    % INTEGRATE_MELT_BY_CATCHMENT
    %
    % melt = integrate_melt_by_catchment(moulin_indices, catchments, node_areas, catchment_melt)

    n_moulins = length(moulin_indices);
    melt = zeros(size(catchments));

    for i=1:length(moulin_indices)
        mask = catchments==(i-1);
        catch_area = node_areas(mask);
        melt(moulin_indices(i)) = sum(catch_area.*catchment_melt(mask));
    end
end