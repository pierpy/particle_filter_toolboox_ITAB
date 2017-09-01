t = 1:600;
x0 = [0 0 1 pi/2];
f = @(t, x) [1; x(3)*sin(x(1)*2*pi*0.03+x(4)); 0; 0];
g = @(t, x) [x(2)];

noise.type='Gaussian';
noise.sigma_w = 0;
noise.sigma_v = 0;
[x,y]=SimSyS(t,x0,f,g,noise);
pf.noise  = noise;
pf.np = 5000;
pf.noise.sigma_w = [0.1 0.1 0.01 0.01];
pf.noise.sigma_v = 0.1;
pf.resampling_strategy = 'multinomial_resampling';
x0 = [0 0 0.8 pi/2-0.1];
[xh,yh,xp,yp,wp]=EstSys(t,x0,y,f,g,pf);
PlotPF_Res(t,x,xh,xp,wp,y,yh)