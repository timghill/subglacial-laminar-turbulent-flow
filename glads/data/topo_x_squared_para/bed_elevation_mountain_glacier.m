function bed = bed_elevation_mountain_glacier(xy, time, pa)
% out = bed_elevation_mountain_glacier(xy, time)
%
% return bed elevation. Maintain time dependence since the GlaDS code
% assumes this is a function of time

dz = 900;
min_thick = 40;

surface = @(x, y) 100*(x+200).^(1/4) + 1/60*x - 2e10^(1/4);

pa_bench = 0.05;
eps = 1e-16;
f = @(x, pa) (surface(6e3, 0) - pa*6e3)/6e3.^2 * x.^2 + pa*x;
g = @(y) 0.5e-6*abs(y).^3;
h = @(x, pa) (-4.5*x/6e3 + 5).*(surface(x,0) - f(x, pa))./(surface(x, 0) - f(x, pa_bench) + eps);

bed_fun = @(x, y, pa) f(x, pa) + g(y).*h(x, pa);

bed = bed_fun(xy(:, 1), xy(:, 2), pa) + dz;

end
