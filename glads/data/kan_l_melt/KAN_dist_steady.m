function melt=KAN_dist_steady(time, pin, dmesh)
    % melt = KAN_dist_steady(time, pin, dmesh)
    % computes steady moulin inputs based on KAN_L AWS temperatures.
    %
    % Uses average melt rate from KAN_mountain_glacier_mean_melt.txt
    %
    % See also KAN_mountain_glacier_compute_average_melt

    ramp = max(0, min(time/86400/365/25 - 1, 1));

    steady_melt = readmatrix('KAN_mountain_glacier_mean_melt.txt');
    % steady_melt = 1.158e-6 + 0*dmesh.tri.nodes(:, 1);

    melt = ramp.*steady_melt + 0.1/86400/365;
end

