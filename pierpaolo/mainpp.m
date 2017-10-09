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

t = (1:60); % simulation time
u = zeros(1,length(t)); % deterministic inplut
u(1,[10 40]) = .1; 
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
noise.sigma_y = [0.001 0.001 0.001 0.001]; % observation noise

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


p = p0;
p(4) = 2.9;
model.noise = noise; 
model.labels = labels;
model.resampling_strategy = 'systematic_resampling';
model.noise.type_x = 'Normal';
model.y = y;
model.xTrue = x;
model.f = f;
model.g = g;
model.T = length(t);                          % signal length
model.Dt = mean(diff(t));               % sampling time
model.nx = length(x0); 
                % number of observation
model.npar = length(p0);                % number of parameters
model.niter = 200;
model.x0 = x0;
model.p0 = p0;
model.f = f;
model.g = g;
model.Q = diag([0.01 0.01 0.01 0.01]);
model.R = diag([0.001 0.001]);
model.u = u;
model.np = 200; % number of particles
model.t = t;
model.stepSize = 0.01;
%[ xh , loglikelihoodProposed ]  = particlefilter( t, x0, p, u, y, f, g, pf );
% [ loglikelihood, theta, proposedThetaAccepted ] = pmh2( t, x0, p, labels, niter, u, y, f, g, pf, 0.01 );
% model = pFilter( t, model );
% model = pSmoother( model );
[ loglikelihood, theta, proposedThetaAccepted ] = pmh2( p, model);
% for i = 1: model.ny
%     figure
%     plot(model.y(i,:))
% end
% 
% 
% for i = 1: model.nx
%     figure
%     plot(model.xTrue(i,:))
%     hold on
%     plot(model.filter.xFilt(i,:))
%     hold on
%     plot(model.smoother.xSmoot(i,:))
% end
