clear 
close all
clc

epsilon = 1; % neuronal efficacy
kas = 0.2632; % signal decay
kaf = 0.1297; %autoregulation
tau0 = 3.14; %transit time
alpha = 0.3663; %stiffness
E0 = 0.4865; % resting Oxygen extraction
%lambda=-0.3; % autocorrelation neuronal

t = (1:120);
u = zeros(1,length(t));
u(1, [5: 30: end]) = 10^-3;
x0 = [0 1 1 1];

p0 = zeros(6,1);
p0(1) = epsilon;
p0(2) = kas;
p0(3) = kaf;
p0(4) = tau0;
p0(5) = alpha;
p0(6) = E0;

f = @balloonInput;
g = @(t, p, x, u) [ x(3); x(4)];
% g = @(t, p, x, u) [ x(1); x(2); x(3); x(4)];

pf.np = 500;
noise.sigma_x = [0 0 0 0 0];
noise.sigma_y = [0.0001 0.0001 0.0001 0.0001];

noise.type_x = 'Normal';
noise.type_y = 'Normal';
noise.type_p = 'Normal';
pf.resampling_strategy = 'systematic_resampling';

[x,y] = SimSys(t, x0, p0, u, f, g, noise);

labels = zeros(6,1);
labels(1) = 0;
labels(2) = 0;
labels(3) = 0;
labels(4) = 1;
labels(5) = 1;
labels(6) = 0;
noise.sigma_p = [ 0.1 0.01 0.01 0.1 0.01 0.01 ];


niter = 100;
% p = p0;
% p(1) = epsilon;
% p(2) = kas;
% p(3) = kaf;
% p(4) = tau0/1.1;
% p(5) = alpha/1.1;
% p(6) = E0;

pf.noise = noise;
pf.noise.type_x = 'Normal';
pf.noise.sigma_x = [ 0.001 0.001 0.001 0.001 ];
noise.sigma_y = [ 0.0001 0.0001 0.0001 0.0001 ];

p4values = 0.1:0.01:0.6;

for p1 = 1: size(p4values,2)
    p1
        p0(2) = p4values(p1);
        [ xh, yh, uwxp, wxp, xp ] = particlefilter( t, x0, p0, u, y, f, g, pf );
        loglikelihood(p1) = sum(log((sum(uwxp,2))));
        clear uwxp

end


