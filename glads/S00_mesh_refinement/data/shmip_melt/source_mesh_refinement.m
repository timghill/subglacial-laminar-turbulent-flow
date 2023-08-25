function melt = source_moulin_shmip_steady(time, pin, dmesh)
    % melt = source_moulin_shmip_steady(time, pin, dmesh, ii_moulin, catchmap)
    % compute steady moulin inputs using SHMIP melt parameterization
    %
    % Uses 25 year winter steady-state + 25 year linear ramp-up of
    % melt intensity to ensure model stability.
    %
    % Uses average melt rate from compute_average_melt
    %
    % See also compute_average_melt

    ramp = max(0, min(time/86400/365/25 - 1, 1));

    %% Load the original dmesh
    ref_meshes = load('../data/mesh/mesh.mat');
    ref_dmesh = ref_meshes.meshes{4};

    n_moulin = 68;
    ref_moulindata = readmatrix(sprintf('./data/shmip_melt/moulins_%03d.txt', n_moulin));
    catchmap = readmatrix(sprintf('./data/shmip_melt/catchment_map_%03d.txt', n_moulin));
    ii_moulin = ref_moulindata(:, 1) + 1;

    ref_xy = ref_dmesh.tri.nodes;

    % Read steady surface melt
    ref_steady_melt = readmatrix('SHMIP_adj_mean_melt.txt');

    ref_area = ref_dmesh.tri.area_nodes;
    ref_catch_melt = integrate_melt_by_catchment(ii_moulin, catchmap, ref_area, ref_steady_melt);

    % Add this to nodes computed for the new mesh
    moulindata = readmatrix(sprintf('../data/moulins/mesh_refinement/moulins_%05d.txt', dmesh.tri.n_nodes));
    moulin_index = moulindata(:, 1) + 1;
    catch_melt = zeros(dmesh.tri.n_nodes, 1);
    for jj=1:n_moulin
        m_in = ref_catch_melt(ii_moulin(jj));
        catch_melt(moulin_index(jj)) = catch_melt(moulin_index(jj)) + m_in;
    
    melt = catch_melt.*ramp;
    % % Load the original dmesh
    % ref_meshes = load('../data/mesh/mesh.mat');
    % ref_dmesh = ref_meshes.meshes{4};

    % n_moulin = 68;
    % moulindata = readmatrix(sprintf('./data/shmip_melt/moulins_%03d.txt', n_moulin));
    % catchmap = readmatrix(sprintf('./data/shmip_melt/catchment_map_%03d.txt', n_moulin));
    % ii_moulin = moulindata(:, 1) + 1;

    % moulins = zeros(ref_dmesh.tri.n_nodes, 1);
    % moulins(ii_moulin) = 1;

    % ref_xy = ref_dmesh.tri.nodes;

    % % Read steady surface melt
    % ref_steady_melt = readmatrix('SHMIP_adj_mean_melt.txt');

    % ref_area = ref_dmesh.tri.area_nodes;
    % ref_catch_melt = integrate_melt_by_catchment(ii_moulin, catchmap, ref_area, ref_steady_melt);

    % % Now find the nearest nodes to inject melt
    % ii_interp = [];
    % xy = dmesh.tri.nodes;
    % catch_melt = zeros(size(xy, 1), 1);
    % for ii=1:length(ref_catch_melt)
    %     meltii = ref_catch_melt(ii);
    %     if meltii>0
    %         dist = sqrt((ref_xy(ii, 1) - xy(:, 1)).^2 + (ref_xy(ii, 2) - xy(:, 2)).^2);
    %         dist(catch_melt>0) = 1e10;
    %         [mindist, ii_min] = min(dist);
    %         catch_melt(ii_min) = meltii;
	%     ii_interp = [ii_interp, ii_min];
    %     end
    % end
    % n_nodes = length(dmesh.tri.nodes);
    % writematrix(ii_interp, sprintf('moulins_%05d.txt', n_nodes))
    % melt = catch_melt.*ramp;
end
