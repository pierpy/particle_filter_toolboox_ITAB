clear 
close all
clc
t = (1:60);
x0 = [0 1 1 1 0];
p0 = [];
f = @balloon;
g = @(t,p, x,u) [x(3); x(4); x(5)];


pf.np = 200;
noise.sigma_x = [0 0 0 0 0];
noise.sigma_y = [ 0.01 0.01 0.001];
noise.sigma_p = [];
noise.type_x = 'Normal';
noise.type_y = 'Normal';
noise.type_p = 'Normal';
pf.resampling_strategy = 'systematic_resampling';

[x,y] = simulatesystem(t, x0, p0, [], f, g, noise);

p = [];
u = [];
pf.noise = noise;
pf.noise.type_x = 'Normal';
pf.noise.sigma_x = [0.01 0.01 0.01 0.01 0.001];
[xh, yh, ph, pp, xp, yp, wpp, wxp] = estimatesystem(t, x0, p, u, y, f, g, pf);
PlotPF_Res(t, p0, ph, pp, wpp, x, xh, xp, wxp, y, yh)
