% This script contains the code necessary for running the nonlinear example 
% in Section 10.2 of the paper
%
% Thomas B. Sch?n, Adrian Wills and Brett Ninness. System Identification 
% of Nonlinear State-Space Models. Automatica, 47(1):39-49, January 2011. 
%

clear;
%====================================================
%===   Simulate the nonlinear state-space model   ===
%====================================================
            % Number of data
m.nstates = 4;
m.nobs = 2; 
m.ninputs = 1; 
tt = (1:600)/10;
N    = length(tt);
u = zeros(1,length(tt));
u(1,[100 400]) = 0.1;
x0 = [0 1 1 1];
m.Q = diag([0 0 0 0]);
m.R = diag([0.000001 0.000001]);
m.P1 = 1;  % Initial state (fully known)
m.X1 = x0';
x = zeros(m.nstates, N+1); 
y = zeros(m.nobs, N);
x(:, 1) = x0';
v = m.Q*randn(m.nstates, N);    % Process noise sequence
e = m.R*randn(m.nobs, N);% Measurement noise sequence

epsilon=1; % neuronal efficacy
kas=0.2632; % signal decay
kaf=0.1297; %autoregulation
tau0=1; %transit time
alpha=0.3663; %stiffness
E0=0.4865; % resting Oxygen extraction

theta = [epsilon kas kaf tau0 alpha E0];
f = @balloonInput;
g = @(t, p, x, u) [ x(3); x(4) ];
m.dt = mean(diff(tt));
for t = 1 : N
    
  x(:, t+1) = x(:, t) + (f( t, theta, x(:, t), u(:, t) ) + v(:, t))*m.dt;
  y(:, t) = g( t, theta, x(:, t), u(:, t) ) + e(:, t);
  
end

% figure
% plot(x(1,:))
% figure
% plot(x(2,:))
% hold on
% plot(x(3,:))
% hold on
% plot(x(4,:))
% figure
% plot(y(1, :))
% hold on 
% plot(y(2, :))

particles = 200;
z.y     = y;
z.xTrue = x(:, 1:N);
m.Q = diag([0.01 0.01 0.01 0.01]); % aggiungo rumore sul modello, non l' ho fatto in simulazione

gPF = ppff(m, particles, z, u, f, g, theta); % Particle filter
% 
% gPS = ppss(m, 200, gPF, u, f, g, theta);    % Particle smoother
% % figure
% % plot(z.y(1,:))
% % hold on
% % plot(z.y(2,:))
% figure
% plot(z.xTrue(1,:))
% hold on
% plot(gPF.Xf(1, :),'g')
% hold on
% plot(gPS.Xs(1, :),'r')
% 
% figure
% plot(z.xTrue(2,:))
% hold on
% plot(gPF.Xf(2,:),'g')
% hold on
% plot(gPS.Xs(2, :),'r')
% 
% figure
% plot(z.xTrue(3,:))
% hold on
% plot(gPF.Xf(3,:),'g')
% hold on
% plot(gPS.Xs(3, :),'r')
% 
% figure
% plot(z.xTrue(4,:))
% hold on
% plot(gPF.Xf(4,:),'g')
% hold on
% plot(gPS.Xs(4, :),'r')

% figure
% plot(gPS.Xs(1, :), 'LineWidth',2)
% legend('measured', 'true', 'filtered', 'smoothed')

%==================================================
%===   EM algorithm using a particle smoother   ===
%==================================================
% opt.miter = 50;     % Number of EM iterations
% 
% % Initialize randomly within 50% of the true value for all parameters
% % a0 = -0.1  + (0.1*0.1)*(-1 + 100*rand(1));
% 
% theta0 = [epsilon kas kaf tau0 alpha E0];    % Initial guess for the parameters
% % Note that we are estimating the square root factor of the covariance in
% % order to enforce the fact that the covariance should be positive
% % semi-definite.
% disp(['Initial parameter guess, thta0: ' num2str(theta0)])
% M = 200;             % Number of particles
% mEst = m;
% 
% mEst.X1 = m.X1;            % Initial mean values known
% m.P1 = 1;               % Spread the initial particles a bit
% theta = theta0;
% theta_est(1, :) = theta0; 
% %options      = optimset('Display','off','LargeScale','off');
% options      = optimset('LargeScale','off', 'Display','on','TolFun',1e-5,'TolX',1e-5); 
% for k = 2: opt.miter
%   % E step
%   gPF = ppff(m, M, z, u, f, g, theta); % Particle filter
%   gPS = ppss(m, M, gPF, u, f, g, theta);    % Particle smoother
%   % M step
%   [theta, Qval] = fminunc('computeQ', theta, options, gPF.xPResampled, z.y, m.Q, m.R, u, f, g, m.dt);
%   theta_est(k, :) = theta;  
%   likelihoodS(k) = sum(sum(gPS.w));
%   likelihoodF(k) = sum(sum(gPF.w));
%   disp(['Iteration nr: ' num2str(k) ' Estimates, theta: ' num2str(theta)])
% end;
