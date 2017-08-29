%% Particle Filter State Space Model Evaluation
% t= time vector
% x0=inizial state
% y= observations
% f=state transition function
% g= obesrvation function
%pf= particle filter structure
% pf.np= number of particles
% pf.noise.type='Gaussian';
%pf.noise.sigma_w= std of state gaussian noise


function [xh,yh,xp,yp,wp]=SimEstimate(t,x0,y,f,g,pf)



T=length(t);                    % signa length
Dt=mean(diff(t));               % sampling time
nx=length(x0);                  % number of states
y0=g(t(1),x0);               % inizial observation
ny=length(y0);                  % number of observation

np = pf.np;                              % number of particles
% initialization of state, observation and noise
xh=zeros(nx,T);
yh=zeros(ny,T);
xp=zeros(nx,T,np);
yp=zeros(ny,T,np);
wp=zeros(T,np);


% Noise generator function
switch pf.noise.type
    
    case 'Gaussian'
        gen_sys_noise = @(z) normrnd(0,z);
        gen_obs_noise = @(z) normrnd(0,z);
        p_obs_noise   = @(v,u,z) normpdf(v,u,z);
        
        
end

% initial state estimate

xh(:,1)=x0;
yh(:,1)=y0;


for i=1:np                          % simulate initial particles
    xp(:,1,i)=x0+gen_sys_noise(pf.noise.sigma_w);
end
wp(:,1,:) = repmat(1/np, np, 1);           % all particles have the same weight

   
wk=1/np*ones(np,1);   
%% Estimate state
for k = 2:T
   fprintf('Time point = %d/%d\n',k,T);
   
   
   for i=1:np
   xp(:,k,i)=xp(:,k-1,i)+f(k,xp(:,k-1,i))*Dt;
   yp(:,k,i)=g(k,xp(:,k,i));
   % weights (when using the PRIOR pdf): eq 63, Ref 1
   wp(k,i) = wk(i) * p_obs_noise (yp(:,k,i),y(:,k),pf.noise.sigma_v);  %%%% mean of p
   end


%% Normalize weight vector
wk=wp(k,:)./sum(wp(k,:),2);
xk=squeeze(xp(:,k,:));

%% Resampling
% remove this condition and sample on each iteration:
[xp(:,k,:), wk] = resample(xk, wk,pf.resampling_strategy);
 disp('Resampling ...')

%% Compute estimated state

   xh(:,k,:) = xh(:,k-1,:) + wk*xk(:,k,:);

%% Store new weights and particles
pf.w(:,k) = wk;
pf.particles(:,:,k) = xk;



xk = xk(:,idx);                    % extract new particles
wk = repmat(1/Ns, 1, Ns);          % now all particles have the same weight
   
   
   
   
   
   
   
   % filtered observation
   yh(:,k) = obs(k, xh(:,k), 0);
end





%% Resampling function
function [xk, wk, idx] = resample(xk, wk, resampling_strategy)

Ns = length(wk);  % Ns = number of particles

% wk = wk./sum(wk); % normalize weight vector (already done)

switch resampling_strategy
   case 'multinomial_resampling'
      with_replacement = true;
      idx = randsample(1:Ns, Ns, with_replacement, wk);
%{
      THIS IS EQUIVALENT TO:
      edges = min([0 cumsum(wk)'],1); % protect against accumulated round-off
      edges(end) = 1;                 % get the upper edge exact
      % this works like the inverse of the empirical distribution and returns
      % the interval where the sample is to be found
      [~, idx] = histc(sort(rand(Ns,1)), edges);
%}
   case 'systematic_resampling'
      % this is performing latin hypercube sampling on wk
      edges = min([0 cumsum(wk)'],1); % protect against accumulated round-off
      edges(end) = 1;                 % get the upper edge exact
      u1 = rand/Ns;
      % this works like the inverse of the empirical distribution and returns
      % the interval where the sample is to be found
      [~, idx] = histc(u1:1/Ns:1, edges);
   % case 'regularized_pf'      TO BE IMPLEMENTED
   % case 'stratified_sampling' TO BE IMPLEMENTED
   % case 'residual_sampling'   TO BE IMPLEMENTED
   case 'percentile_resampling'
      % this is performing latin hypercube sampling on wk
       p=10;
       good_p(:,1)=find(wk>prctile(wk,p));
      % this works like the inverse of the empirical distribution and returns
      % the interval where the sample is to be found

   % case 'regularized_pf'      TO BE IMPLEMENTED
   % case 'stratified_sampling' TO BE IMPLEMENTED
   % case 'residual_sampling'   TO BE IMPLEMENTED
   
   otherwise
      error('Resampling strategy not implemented')
end