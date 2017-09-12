clear 
close all
clc



epsilon=1; % neuronal efficacy
kas=0.2632; % signal decay
kaf=0.1297; %autoregulation
tau0=3.14; %transit time
alpha=0.3663; %stiffness
E0=0.4865; % resting Oxygen extraction
lambda=-0.3; % autocorrelation neuronal






t = (1:120);
x0 = [0 1 1 1 0];

p0=zeros(7,2);
p0(1,1)=epsilon;
p0(2,1)=kas;
p0(3,1)=kaf;
p0(4,1)=tau0;
p0(5,1)=alpha;
p0(6,1)=E0;
p0(7,1)=lambda;
p0(4,2)=1;

f = @balloon;
g = @(t,p, x,u) [x(3); x(4); x(5)];


pf.np = 200;
noise.sigma_x = [0 0 0 0 10^-2];
noise.sigma_y = [0.01 0.01 0.001];
noise.sigma_p = [0 0 0 10^-1 0 0 0];
noise.type_x = 'Normal';
noise.type_y = 'Normal';
noise.type_p = 'Normal';
pf.resampling_strategy = 'systematic_resampling';

[x,y] = SimSys(t, x0, p0, [], f, g, noise);

p= p0;
p(4,1)=tau0/1.5;

u = [];
pf.noise = noise;
pf.noise.type_x = 'Normal';
pf.noise.sigma_x = [0.01 0.01 0.01 0.01 0.01];
noise.sigma_y = [0.01 0.01 0.0001];
[xh, yh, ph, pp, xp, yp, wpp, wxp] = EstSys(t, x0, p, u, y, f, g, pf);
PlotPF_Res(t, p0, ph, pp, wpp, x, xh, xp, wxp, y, yh)
