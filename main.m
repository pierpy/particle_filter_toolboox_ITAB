clear
t = (1:50)*4;
x0 = [0 0];
p0 = [0.05];
% f = @(t,p, x,u) [1; p(1)*sin(x(1)*2*pi*p(2))];
u = zeros(1,200);
u(1,5) = 1;
f = @(t,p, x,u) [1; -p(1)*x(2) + u(1)];
g = @(t,p, x,u) [x(2)];

noise.sigma_x = [0 0];
noise.sigma_y = [0.01];
noise.type_x = 'Normal';
noise.type_y = 'Normal';
noise.type_p = 'Normal';
pf.resampling_strategy = 'systematic_resampling';
[x,y] = SimSys(t, x0, p0, u, f, g, noise);
% figure
% plot(x(2,:))
% figure
% plot(y)
tic
niter = 10000;
pf.noise.type_x = 'Normal';
pf.noise.sigma_x = [0.01 0.01 0.01 0.01];
noise.sigma_y = [0.01];
pf.noise  = noise;
pf.np = 500;
pf.noise.sigma_x = [0.1 0.1];
pf.noise.sigma_y = 0.1;
pf.noise.sigma_p = [0.01];
pf.noise.type='Gaussian';
pf.resampling_strategy = 'multinomial_resampling';
x0 = [0 0];
p=[0.1];

labels = 1;
tic 
[xh, yh, ph, xp, yp, uwxp, wxp, p, p_theta] = EstSys3 (t, x0, p, labels, niter, u, y, f, g, pf);
toc
% PlotPF_Res(t,p0,ph,pp,wpp,x,xh,xp,wxp,y,yh)