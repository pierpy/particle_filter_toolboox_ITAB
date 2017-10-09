%% State Space Model simulation
% t= time vector
% x0 = inizial state
% f= state transition function
% g = obesrvation function
% noise = noise structure
% noise.type = 'Gaussian';
% noise.sigma_w = std of state gaussian noise
% noise.sigma_v = std of observation gaussian noise
% Example:
% 2 states: state 1 equal to time, state 2 equal to sin at 0.03 Hz
% obervation is the state 2
% t=1:60
% x0=[0 0];
% p0= [parameters]
% f = @(t,p, x) [1; p(1)*sin(x(1)*2*pi*0.03)]
% g=@(t,p,x)[x(2)]
function [x,y] = SimSys(t, x0, p0, u, f, g, noise)

T = length(t);                    % simulation length
Dt = mean(diff(t));               % sampling time
nx = length(x0);                  % number of states
if isempty(u)
    y0 = g(t(1), p0(:,1), x0, []);           % initial observation
else
    y0 = g(t(1), p0(:,1), x0, u(:,1));    % initial observation
end              % inizial observation
ny = length(y0);                  % number of observation
% initialization of state, observation and noise
x = zeros(nx,T);               
y = zeros(ny,T);
w = zeros(nx,T);  
v = zeros(ny,T);
x(:,1) = x0;
y(:,1) = y0;

% Simulate system
for k = 2:T
    for i = 1: nx
        w(i,k) = random (noise.type_x, 0, noise.sigma_x(i)); 
    end% simulate process noise
    for i = 1: ny
        v(i,k) = random (noise.type_y, 0, noise.sigma_y(i));    
    end% simulate observation noise
   if isempty(u)
       x(:,k) = x(:,k-1) + (f(k,p0(:,1), x(:,k-1),[]) + w(:,k))*Dt;     % simulate state
       y(:,k) = g(k,p0(:,1), x(:,k),[])+v(:,k);     % simulate observation
   else
       x(:,k) = x(:,k-1) + (f(k, p0(:,1), x(:,k-1), u(:,k)) + w(:,k))*Dt;
       y(:,k) = g(k,p0(:,1), x(:,k), u(:,k)) +v(:,k);     % simulate observation
   end
end