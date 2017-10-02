clear
t = (1:100);
x0 = [0];
p0 = [-.1]';
u = zeros(1,100);
u(1,10) = 1;
% f = @(t, p, x, u) [1; p(1)*sin(2*pi*p(2)*x(1))];
f = @(t, p, x, u) [p(1)*x(1) + u(1)];
% u = zeros(1,200);
% u(1,5) = 1;
% f = @(t,p, x,u) [1; -p(1)*x(2) + u(1)];
g = @(t,p, x, u) [ x(1) ];
noise.sigma_x = [0];
noise.sigma_y = [.1];
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
pf.np = 500;
pf.noise.sigma_x = [0.1];
pf.noise.sigma_y = .1;
pf.noise.sigma_p = [10^-2 1 10^-1 10^-3];
pf.noise.type = 'Gaussian';
pf.resampling_strategy = 'multinomial_resampling';
x0 = [0];


% tic 
% p1values = -0.00001:-0.01:-0.5;

% p2values = 0.01:0.01:0.1;
% for p1 = 1:size(p1values,2)
%     p1
%     for p2 = 1:size(p2values,2)
%         p2
%         p0(1) = p1values(p1);
%         p0(5) = p2values(p2);
%         [ xh, yh, uwxp, wxp, xp ] = particlefilter( t, x0, p0, u, y, f, g, pf );
%         [ xhs, yhs, xps, wps, uwps ] = backwardsmoother( t, p0, u, f, wxp, xp, xh, yh, pf );
%         uwps = squeeze(uwps);
%         wps = squeeze(wps);
%         loglikelihood(p1) = sum(log((sum(uwxp,2))));
%         clear uwxp
%     end
% end
[ xh, yh, uwxp, wxp, xp, xpnr ] = particlefilter( t, x0, p0, u, y, f, g, pf );
[ J, xps, xhs] = particlesmoother( t, p0, u, f, xp, pf );
%  if sum(labels > 0)
%      [ loglikelihood, ph ] = pmh( t, x0, p, labels, niter, u, y, f, g, pf );
%  end
% toc
% % PlotPF_Res(t,p0,ph,pp,wpp,x,xh,xp,wxp,y,yh)
figure
plot(x)
hold on
plot(xh, 'r')
hold on
plot(xhs, 'g')

rf = rms(x-xh)
rs = rms(x-xhs)