











t = (1:2000)/10;
x0 = [0 1 1 1 0];
f = @balloon;
g = @(t,p, x,u) [x(3); x(4); x(5)];
noise.type='Gaussian';
noise.sigma_x = 0;
noise.sigma_y = 10^-2;
[x,y]=SimSyS(t,x0,1,[],f,g,noise);
pf.noise  = noise;
pf.np = 500;
pf.noise.sigma_x=[10^-2 10^-2 10^-2 10^-2 10^-2];
pf.noise.sigma_y=10^-2;
pf.noise.sigma_p=[];
pf.noise.type='Gaussian';
pf.resampling_strategy = 'multinomial_resampling';
x0 = [0 1 1 1 0];
p=[];
u=[];
[xh,yh,ph,pp,xp,yp,wpp,wxp]=EstSys(t,x0,p,u,y,f,g,pf);
PlotPF_Res(t,p0,ph,pp,wpp,x,xh,xp,wxp,y,yh)


