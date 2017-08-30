t = 1:600;
x0 = [0 0];
f = @(t, x) [1; sin(x(1)*2*pi*0.03)];
g = @(t, x) [x(2)];

noise.type='Gaussian';
noise.sigma_w = 0;
noise.sigma_v = 3;
%[x,y]=SimSyS(t,x0,f,g,noise);
pf.noise  = noise;
pf.np = 500;
pf.noise.sigma_w = [0.1 0.1];
pf.noise.sigma_v = 1000;
pf.resampling_strategy = 'multinomial_resampling';

[xh,yh,xp,yp,wp]=EstSys(t,x0,y,f,g,pf);
PlotPF_Res(t,x,xh,xp,wp,y,yh)