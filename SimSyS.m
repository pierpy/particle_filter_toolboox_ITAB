%% State Space Model simulation
% t= time vector
% x0=inizial state
% f= state transition function
% g= obesrvation function
% noise= noise structure
% noise.type='Gaussian';
%noise.sigma_w= std of state gaussian noise
%noise.sigma_v==std of observation gaussian noise


% Example:
% 2 states: state 1 equal to time, state 2 equal to sin at 0.03 Hz
% obervation is the state 2
% t=1:60
% x0=[0 0];
% p0= [parameters]
% f = @(t,p, x) [1; p(1)*sin(x(1)*2*pi*0.03)]
% g=@(t,p,x)[x(2)]

function [x,y]=SimSyS(t,x0,p0,f,g,noise)

T=length(t);                    % simulation length
Dt=mean(diff(t));               % sampling time
nx=length(x0);                  % number of states
y0=g(t(1),p0,x0);               % inizial observation
ny=length(y0);                  % number of observation

% initialization of state, observation and noise
x = zeros(nx,T);               
y = zeros(ny,T);
w = zeros(nx,T);  
v = zeros(ny,T);
x(:,1)=x0;
y(:,1)=y0;

% Noise generator function
  switch noise.type
      
      case 'Gaussian'
      gen_state_noise = @(z) normrnd(0, z);   
      gen_obs_noise = @(z) normrnd(0, z);
      
  end
  
% Simulate system
for k = 2:T
 
   w(:,k) = gen_state_noise(noise.sigma_w);              % simulate process noise
   v(:,k) = gen_obs_noise(noise.sigma_v);               % simulate observation noise
   x(:,k) = (x(:,k-1)+f(k,p0, x(:,k-1))+w(:,k))*Dt;     % simulate state
   y(:,k) = g(k,p0, x(:,k))+v(:,k);     % simulate observation
   
end