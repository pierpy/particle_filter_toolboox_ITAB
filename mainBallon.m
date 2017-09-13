clear 
close all
clc



epsilon=1; % neuronal efficacy
kas=0.2632; % signal decay
kaf=0.1297; %autoregulation
tau0=3.14; %transit time
alpha=0.3663; %stiffness
E0=0.4865; % resting Oxygen extraction
%lambda=-0.3; % autocorrelation neuronal






t = (1:600)/10;
u=zeros(1,length(t));
u(1,[100 400])=1;

x0 = [0 1 1 1];

p0=zeros(7,2);
p0(1,1)=epsilon;
p0(2,1)=kas;
p0(3,1)=kaf;
p0(4,1)=tau0;
p0(5,1)=alpha;
p0(6,1)=E0;







f=@balloonInput;
g = @(t,p, x,u) [x(3); x(4); u];



pf.np =300;
noise.sigma_x = [0 0 0 0 0];
noise.sigma_y = [0.01 0.01 0.001];

noise.type_x = 'Normal';
noise.type_y = 'Normal';
noise.type_p = 'Normal';
pf.resampling_strategy = 'systematic_resampling';

[x,y] = SimSys(t, x0, p0, u, f, g, noise);

p= p0;
p(4,1)=tau0/1.2;
%p(5,1)=alpha/1.5;
%p0(2,1)=kas/1.5;
labels=zeros(6,1);
labels(4)=0;
noise.sigma_p = [0 0 0 10^-1 0 0];

tic
niter=10;
pf.noise = noise;
pf.noise.type_x = 'Normal';
pf.noise.sigma_x = [0.01 0.01 0.01 0.01];
noise.sigma_y = [0.01 0.01 0.001];
[xh, yh, ph, xp, yp, uwxp, wxp,p,p_theta] = EstSys3 (t, x0, p,labels,niter, u, y, f, g, pf);
toc

PlotPF_Res(t,p,ph,x,xh,xp,wxp,y,yh)

