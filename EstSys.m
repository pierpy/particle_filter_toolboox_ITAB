%% Particle Filter State Space Model Evaluation
% t= time vector
% x0 = inizial state
% po = inizial parameters first column, estimate or not estimate [1 0]
% second column
% y = observations
% f = state transition function
% g = obesrvation function
% pf = particle filter structure
% pf.np = number of particles
% pf.pf.noise.type ='Gaussian';
% pf.pf.noise.sigma_w= std of state gaussian pf.noise
function [xh, yh, ph, pp, xp, yp, wpp, wxp] = EstSys (t, x0, p0, u, y, f, g, pf)
    T = length(t);                    % signal length
    Dt = mean(diff(t));               % sampling time
    nx = length(x0);                  % number of states
    if isempty(u)
     y0 = g(t(1), p0(:,1), x0, []);           % initial observation
    else  
     y0 = g(t(1), p0(:,1), x0, u(:,1));    % initial observation    
    end
    ny = length(y0);                  % number of observation
    npar = length(p0(:,1));                % number of parameters
    np = pf.np;                              % number of particles
    % initialization of state, observation and pf.noise
    xh = zeros(nx, T);
    yh = zeros(ny, T);
    ph = zeros(npar, T);
    xp = zeros(nx, T, np);
    yp = zeros(ny, T, np);
    pp = zeros(npar, T, np);
    wxp = zeros(T, np);
    wpp = [];
    if ~ isempty(p0)
    % initial parameters estimate
        wpp = zeros(T,np);
        ph(:,1) = p0(:,1);   
        npar = size(ph,1);
        for j =1: npar
            for i = 1: np    % simulate initial particles
                if p0(j,2)==1
                pp(j, 1, i) = p0(j,1) + random (pf.noise.type_p, 0, pf.noise.sigma_p(j));  
                else
                pp(j, 1, i) = p0(j,1);    
                end
            end
        end
        wpp(1,:) = repmat(1/np, np, 1);           % all particles have the same weight
        wpk = 1/np*ones(np, 1);  
    end
    % initial states estimate
    xh(:, 1) = x0;
    yh(:, 1) = y0;
    ph(:, 1) = p0(:, 1);
    for j = 1: nx
        for i = 1: np                          % simulate initial particles
            xp(j,1,i) = x0(j) + random (pf.noise.type_x, 0, pf.noise.sigma_x(j));
        end
    end
    wxp(1,:) = repmat(1/np, np, 1);           % all particles have the same weight
    wxk = 1/np*ones(np, 1);
    % Estimate parameters and states
    for k = 2: T
        fprintf('Time point = %d/%d\n',k,T);
        if ~ isempty(p0)
            % parameters estimates with metropolitan marginal hasting
            for i = 1: size(pp,3)
                w = zeros(npar,1);
                for j=1: npar
                    if p0(j,2)==1
                        w(j) = random (pf.noise.type_p, 0, pf.noise.sigma_p(j));
                    else
                        w(j) =0;  
                    end
                end
                % marginal metropolitan hasting
                pProposed = pp(:,k-1,i) + w;
                if isempty(u)
                    xpold = xp(:, k-1, i) + (f(k, pp(:,k-1,i), xp(:,k-1,i), [])) * Dt;
                    xpnew = xp(:, k-1, i) + (f(k, pProposed, xp(:,k-1,i), [])) * Dt;
                    yold = g(k, pp(:, k-1, i), xpold, []);
                    ynew = g(k, pProposed, xpnew, []);
                else
                    xpold = xp(:, k-1,i) + (f(k, pp(:,k-1,i), xp(:,k-1,i), u(:,k))) * Dt;
                    xpnew = xp(:, k-1,i) + (f(k, pProposed, xp(:,k-1,i), u(:,k))) * Dt;
                    yold = g(k, pp(:,k-1,i), xpold, u(:,k));
                    ynew = g(k, pProposed, xpnew, u(:,k));
                end
                if evalpdfobs( y(:,k), ynew, pf.noise ) > evalpdfobs( y(:,k), yold, pf.noise )
                    wpp(k,i) = wpk(i) * evalpdfobs( y(:,k), ynew, pf.noise );  %%%% mean of p's yp given xp
                    pp(:,k,i) = pProposed;
                else
                    wpp(k,i) = wpk(i) * evalpdfobs( y(:,k), yold, pf.noise );  %%%% mean of p's yp given xp
                    pp(:,k,i) = pp(:,k-1,i);
                end
            end
            % Normalize weight vector
            wpp(k,:) = wpp(k,:)./sum(wpp(k,:), 2);
            wpk = wpp(k,:);
            if npar == 1
                pk=squeeze(pp(:,k,:))';
            else
                pk=squeeze(pp(:,k,:));        
            end
            Neff = 1/sum(wpk.^2);
            resample_percentage = 0.5;
            Nt = resample_percentage * np;
            if Neff < Nt
            % Resampling
                disp('Parameters Particles Resampling ...')
                % remove this condition and sample on each iteration:
                [pp(:,k,:), wpk,~] = resample(pk, wpk', pf.resampling_strategy);
            end

            if size(pp,1) == 1
                ph(:,k) = squeeze(pp(:,k,:))'*wpk';
            else
                ph(:,k) = squeeze(pp(:,k,:))*wpk';
            end
        end
        % states estimates
        for i=1:np
            w=zeros(nx,1);
            for j=1:nx
                w(j) = random (pf.noise.type_x, 0, pf.noise.sigma_x(j));
            end
            if isempty(u)
                xp(:,k,i) = xp(:,k-1,i) + (f(k,pp(:,k,i), xp(:,k-1,i), []) + w) * Dt;
                yp(:,k,i) = g(k,pp(:,k,i), xp(:,k,i), []);
            else
                xp(:,k,i) = xp(:,k-1,i) + (f(k,pp(:,k,i), xp(:,k-1,i), u(:,k)) + w) * Dt;
                yp(:,k,i) = g(k,pp(:,k,i), xp(:,k,i), u(:,k));
            end
            % weights (when using the PRIOR pdf): eq 63, Ref 1
            wxp(k,i) = wxk(i) * evalpdfobs( y(:,k), yp(:,k,i), pf.noise );  %%%% mean of p's yp given xp
            
        end
        % Normalize weight vector
        wxp(k,:) = wxp(k,:)./sum(wxp(k,:),2);
        wxk = wxp(k,:);
        xk = squeeze(xp(:,k,:));        
        Neff = 1/sum(wxk.^2);
        resample_percentage = 0.5;
        Nt  = resample_percentage*np;
        if Neff < Nt
        % Resampling
        disp('States Particles Resampling ...')
        % remove this condition and sample on each iteration:
        [xp(:,k,:), wxk,~] = resample(xk, wxk', pf.resampling_strategy);
       end
        % Compute estimated state
        xh(:,k) = squeeze(xp(:,k,:))*wxk';
        if isempty(u)
            % filtered observation
            yh(:,k) = g(k, ph(:,k),xh(:,k),[]);
        else
            % filtered observation
            yh(:,k) = g(k, ph(:,k),xh(:,k),u(:,k));
        end
    end
end

%% Resampling function
function [xk, wk, idx] = resample(xk, wk, resampling_strategy)
    Ns = length(wk);  % Ns = number of particles
    switch resampling_strategy
        case 'multinomial_resampling'
            with_replacement = true;
            idx = randsample(1:Ns, Ns, with_replacement, wk);
        case 'systematic_resampling'
            % this is performing latin hypercube sampling on wk
            edges = min([0 cumsum(wk)'],1); % protect against accumulated round-off
            edges(end) = 1;                 % get the upper edge exact
            u1 = rand/Ns;
            [~, idx] = histc(u1:1/Ns:1, edges);
        case 'percentile_resampling'
            % this is performing latin hypercube sampling on wk
            p=10;
            good_p(:,1)=find(wk>prctile(wk,p));
        otherwise
            error('Resampling strategy not implemented')
    end
    xk = xk(:,idx);
    wk = repmat(1/Ns,1,Ns);
end