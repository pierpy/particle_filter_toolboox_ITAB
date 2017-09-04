%% Particle Filter State Space Model Evaluation
% t= time vector
% x0=inizial state
% po=inizial parameters
% y= observations
% f=state transition function
% g= obesrvation function
%pf= particle filter structure
% pf.np= number of particles
% pf.noise.type='Gaussian';
%pf.noise.sigma_w= std of state gaussian noise


function [xh,yh,ph,pp,xp,yp,wpp,wxp]=EstSys(t,x0,p0,y,f,g,pf)

    T=length(t);                    % signal length
    Dt=mean(diff(t));               % sampling time
    nx=length(x0);                  % number of states
    y0=g(t(1),p0,x0);               % inizial observation
    ny=length(y0);                  % number of observation
    npar=length(p0);                % number of parameters
    np = pf.np;                              % number of particles
    % initialization of state, observation and noise
    xh=zeros(nx,T);
    yh=zeros(ny,T);
    ph=zeros(npar,T);
    xp=zeros(nx,T,np);
    yp=zeros(ny,T,np);
    pp=zeros(npar,T,np);
    wxp=zeros(T,np);
    wpp=zeros(T,np);


    % Noise generator function
    switch pf.noise.type

        case 'Gaussian'
            gen_state_noise = @(z) normrnd(0,z);
            p_obs_noise   = @(v,u,z) normpdf(v,u,z);
            gen_par_noise = @(z) normrnd(0,z);

    end

    
    % initial parameters estimate
  
    ph(:,1)=p0;
    
    npar=size(ph,1);
    for j=1:npar
        for i=1:np                          % simulate initial particles
            pp(j,1,i)=p0(j)+gen_par_noise(pf.noise.sigma_p(j));           
            
        end
    end
    wpp(1,:) = repmat(1/np, np, 1);           % all particles have the same weight
    wpk=1/np*ones(np,1);
    
    
    % initial states estimate

    xh(:,1)=x0;
    yh(:,1)=y0;
    ph(:,1)=p0;
    
    nx=size(xh,1);
    for j=1:nx
        for i=1:np                          % simulate initial particles
            xp(j,1,i)=x0(j)+gen_state_noise(pf.noise.sigma_w(j));
        end
    end
    wxp(1,:) = repmat(1/np, np, 1);           % all particles have the same weight
    wxk=1/np*ones(np,1);
    
    
    

    %% Estimate parameters and states
    for k = 2:T
        fprintf('Time point = %d/%d\n',k,T);
        
        %% parameters estimates
        for i=1:np
            w=zeros(npar,1);
            for j=1:npar
                w(j)=gen_par_noise(pf.noise.sigma_p(j));
            end
            
            pProposed=pp(:,k-1,i)+w;
            
            xpold=xp(:,k-1,i)+(f(k,pp(:,k-1,i),xp(:,k-1,i)))*Dt;
            xpnew=xp(:,k-1,i)+(f(k,pProposed,xp(:,k-1,i)))*Dt;
            yold=g(k,pp(:,k-1,i),xpold);
            ynew=g(k,pProposed,xpnew);                       
            if mean((p_obs_noise (y(:,k)-ynew,0,pf.noise.sigma_v)+1e-300))>mean((p_obs_noise (y(:,k)-yold,0,pf.noise.sigma_v)+1e-300))
            wpp(k,i) = wpk(i) *mean((p_obs_noise (y(:,k)-ynew,0,pf.noise.sigma_v)+1e-300));  %%%% mean of p's yp given xp
            pp(:,k,i)=pProposed;
            else
            wpp(k,i) = wpk(i) *mean((p_obs_noise (y(:,k)-yold,0,pf.noise.sigma_v)+1e-300));  %%%% mean of p's yp given xp
            pp(:,k,i)=pp(:,k-1,i);
            end
                
            
        end


        % Normalize weight vector
        wpp(k,:)=wpp(k,:)./sum(wpp(k,:),2);
        wpk=wpp(k,:);
        if npar==1
        pk=squeeze(pp(:,k,:))';
        else
        pk=squeeze(pp(:,k,:));        
        end
        
        
     
        Neff=1/sum(wpk.^2);
        resample_percentage=0.5;
        Nt=resample_percentage*np;
        if Neff<Nt
        % Resampling
        disp('Parameters Particles Resampling ...')
        % remove this condition and sample on each iteration:
        [pp(:,k,:), wpk,~] = resample(pk, wpk','systematic_resampling');
%         [~,idx]=sort(wk);
%         idx = randsample(1:np, np, with_replacement, wk(idx(idx/2+1:end)));
        end

        if size(pp,1)==1
        ph(:,k) = squeeze(pp(:,k,:))'*wpk';
        else
        ph(:,k) = squeeze(pp(:,k,:))*wpk';
        end
        
        %% state estimates
        for i=1:np
            w=zeros(nx,1);
            for j=1:nx
                w(j)=gen_state_noise(pf.noise.sigma_w(j));
            end
            xp(:,k,i)=xp(:,k-1,i)+(f(k,pp(:,k,i),xp(:,k-1,i))+w)*Dt;
            yp(:,k,i)=g(k,pp(:,k,i),xp(:,k,i));
            % weights (when using the PRIOR pdf): eq 63, Ref 1
            wxp(k,i) = wxk(i) *mean((p_obs_noise (y(:,k)-yp(:,k,i),0,pf.noise.sigma_v)+1e-300));  %%%% mean of p's yp given xp
            
        end


        % Normalize weight vector
        wxp(k,:)=wxp(k,:)./sum(wxp(k,:),2);
        wxk=wxp(k,:);
        xk=squeeze(xp(:,k,:));

        
        Neff=1/sum(wxk.^2);
        resample_percentage=0.5;
        Nt=resample_percentage*np;
        if Neff<Nt
        %% Resampling
        disp('States Particles Resampling ...')
        % remove this condition and sample on each iteration:
        [xp(:,k,:), wxk,~] = resample(xk, wxk','systematic_resampling');
%         [~,idx]=sort(wk);
%         idx = randsample(1:np, np, with_replacement, wk(idx(idx/2+1:end)));
       end
        %% Compute estimated state

        xh(:,k) = squeeze(xp(:,k,:))*wxk';

        % filtered observation
        yh(:,k) = g(k, ph(:,k),xh(:,k));
    end
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
      xk=xk(:,idx);
       wk=repmat(1/Ns,1,Ns);
end