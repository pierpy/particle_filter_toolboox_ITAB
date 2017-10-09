% attempt PMH for parameter estimation. 
% first the system is simulated. 
% then the call to pmh function perform PMH algorithm:
% pmh function calls  self explamatory particlefilter function

clear 
close all
clc

% true model parameters
epsilon = 1; % neuronal efficacy
kas = 0.2632; % signal decay
kaf = 0.1297; %autoregulation
tau0 = 3.14; %transit time
alpha = 0.3663; %stiffness
E0 = 0.4865; % resting Oxygen extraction
%lambda=-0.3; % autocorrelation neuronal

t = (1:600)/10; % simulation time
u = zeros(1,length(t)); % deterministic inplut
u(1,[100 400]) = 0.1; 
x0 = [0 1 1 1]; % initial state 

p0 = zeros(6,1); % create true parameters vector
p0(1) = epsilon;
p0(2) = kas;
p0(3) = kaf;
p0(4) = tau0;
p0(5) = alpha;
p0(6) = E0;

f = @balloonInput; % state evolution function
g = @(t, p, x, u) [ x(3); x(4) ]; % observation function
% g = @(t, p, x, u) [ x(1); x(2); x(3); x(4)];

% noise structure for simulaton 
noise.sigma_x = [0 0 0 0 0]; % no noise on the simulated states
noise.sigma_y = [0.00001 0.00001 0.00001 0.00001]; % observation noise

% distribution noise type
noise.type_x = 'Normal'; 
noise.type_y = 'Normal';
noise.type_p = 'Normal';
% simulate sistem
[x,y] = SimSys(t, x0, p0, u, f, g, noise);

labels = zeros(6,1); % flag for which parameter to be estimated
labels(1) = 0;
labels(2) = 0;
labels(3) = 0;
labels(4) = 1;
labels(5) = 0;
labels(6) = 0;

noise.sigma_p = [ 0.1 0.01 0.01 0.1 0.01 0.01 ]; % noise variance on parameters


niter = 100;
p = p0;
pf.noise = noise; % structure for particlefilter function
pf.resampling_strategy = 'systematic_resampling';
pf.np = 200; 
pf.noise.type_x = 'Normal';
pf.noise.sigma_x = [ 0.001 0.001 0.001 0.001 ]; % state evolution noise
pf.R = diag([0.000001 0.000001]);
pf.Q = diag([0.001 0.01 0.01 0.001]);
noise.sigma_y = [ 0.00001 0.00001 0.00001 0.00001 ]; % observation noise

%[ xh , loglikelihoodProposed ]  = particlefilter( t, x0, p, u, y, f, g, pf );
% [ loglikelihood, theta, proposedThetaAccepted ] = pmh2( t, x0, p, labels, niter, u, y, f, g, pf, 0.01 );
tic
[ xh , loglikelihoodProposed, particles ]  = pFilter( t, x0, p0, u, y, f, g, pf );
 toc
str = ppss(pf.Q, pf.np, particles, u, f, g, p);
% check if loglikelihood follows parameters.
% p4values = 0.1:0.01:0.6;
% 
% for p1 = 1: size(p4values,2)
%     p1
%         p0(2) = p4values(p1);
%         [ xh, yh, uwxp, wxp, xp ] = particlefilter( t, x0, p0, u, y, f, g, pf );
%         loglikelihood(p1) = sum(log((sum(uwxp,2))));
%         clear uwxp
% 
% end


