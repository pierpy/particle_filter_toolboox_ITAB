clear
t = (1:100);
x0 = [0];
p0 = [0.5 25 8 0.05]';
% u = zeros(1,100);
% u(1,5) = 1;
u = [];
% f = @(t, p, x, u) [1; p(1)*sin(2*pi*p(2)*x(1))];
% f = @(t, p, x, u) [p(1)*x(1) + u(1)];

f = @(t, p, x, u) [(p(1)-1)*x(1) + p(2)*x(1) / (1 + x(1)^2) + p(3)*cos(1.2*t)];
% u = zeros(1,200);
% u(1,5) = 1;
% f = @(t,p, x,u) [1; -p(1)*x(2) + u(1)];

g = @(t,p, x, u) [ p(4)*x(1)^2 ];
noise.sigma_x = [0];
noise.sigma_y = [.01];
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
niter = 200;
pf.noise.type_x = 'Normal';
pf.noise  = noise;
pf.np = 100;
pf.noise.sigma_x = [0.01];
pf.noise.sigma_y = .01;
pf.noise.sigma_p = [10^-3 ];
pf.noise.type = 'Gaussian';
pf.resampling_strategy = 'multinomial_resampling';
x0 = [0];
% 
% [ xh, yh, uwxp, wxp, xp, xpnr ] = particlefilter( t, x0, p0, u, y, f, g, pf );
% [ J, xps, xhs] = particlesmoother( t, p0, u, f, xp, pf );

p0 =[0.4 20 6 0.04]';
options      = optimset('LargeScale','off', 'Display','off','TolFun',1e-5,'TolX',1e-15); 
for k = 1: 10
    
  % E step
  [ xh, yh, uwxp, wxp, xp, xpnr ] = particlefilter( t, x0, p0, u, y, f, g, pf );
  [ J, xps, xhs] = particlesmoother( t, p0, u, f, xp, pf );

  % M step
  [p0, Qval] = fminunc('computeQfunexp', p0, options, xps, y, f, g, u, pf);
 
  mEstStore(:,k) = p0;    % Store the iterates
  disp(['Iteration nr: ' num2str(k) ' Estimates, p0: ' num2str(p0(1))])
end;


% p4values = -0.00001:-0.01:-0.5;
% 
% for p1 = 1: size(p4values,2)
%     p1
%     p0 = p4values(p1);
%     [ xh, yh, uwxp, wxp, xp, xpnr ] = particlefilter( t, x0, p0, u, y, f, g, pf );
%     [ J, xps, xhs] = particlesmoother( t, p0, u, f, xp, pf );
%     Qfun = computeQfunexp(p0, xps, y, f, g, u, pf);
%     QQ(p1) = Qfun;
% end
