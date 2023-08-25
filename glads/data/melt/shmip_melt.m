function melt = shmip_melt(t)
% SHMIP_MELT : load mean melt array
%

fpath = 'SHMIP_mean_melt.txt';

melt = readmatrix(fpath);

end

