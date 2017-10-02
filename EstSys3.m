%% Particle Filter State Space Model Evaluation with 
% t= time vector
% x0 = inizial state
% p0 = inizial parameters 
% second column
% y = observations
% f = state transition function
% g = obesrvation function
% pf = particle filter structure
% pf.np = number of particles
% pf.pf.noise.type ='Gaussian';
% pf.pf.noise.sigma_w= std of state gaussian pf.noise
% niter=number of parameter estimation
% labels = length of p0 vector with 1 estimate 0 do not estimate


function [xh, yh, ph, xp, yp, uwxp, wxp, p, p_theta] = EstSys3 (t, x0, p0, labels, niter, u, y, f, g, pf)

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
    wxp = zeros(T, np);
    uwxp = zeros(T,np);
    
    % initial states estimate
    xh(:, 1) = x0;
    yh(:, 1) = y0;
    for j = 1: nx
        for i = 1: np                          % simulate initial particles
            xp(j,1,i) = x0(j) + random (pf.noise.type_x, 0, pf.noise.sigma_x(j));
        end
    end
    wxp(1,:) = repmat(1/np, np, 1); 
    uwxp(1,:) = repmat(1/np, np, 1); % all particles have the same weight
    wxk = 1/np*ones(np,1);
   
    for k=2:T
        fprintf('Time point = %d/%d\n',k,T);
        % states estimates
        for i=1: np
            w=zeros(nx,1);
            for j=1:nx
                w(j) = random (pf.noise.type_x, 0, pf.noise.sigma_x(j));
            end
            if isempty(u)
                xp(:,k,i) = xp(:,k-1,i) + (f(k, p0, xp(:,k-1,i), []) + w) * Dt;
                yp(:,k,i) = g(k, p0, xp(:,k,i), []);
            else
                xp(:,k,i) = xp(:,k-1,i) + (f(k, p0, xp(:,k-1,i), u(:,k)) + w) * Dt;
                yp(:,k,i) = g(k, p0, xp(:,k,i), u(:,k));
            end
            % weights (when using the PRIOR pdf): eq 63, Ref 1
            wxp(k,i) = wxk(i) * evalpdfobs( y(:,k), yp(:,k,i), pf.noise );  %%%% mean of p's yp given xp
            
        end
        % Normalize weight vector
        uwxp(k,:) = wxp(k,:);
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
            yh(:,k) = g(k, p0, xh(:,k),[]);
        else
            % filtered observation
            yh(:,k) = g(k, p0, xh(:,k), u(:,k));
        end
    end

    p_theta1 = sum(log((sum(uwxp,2))));
    p = [];
    p_theta = [];
    if sum(labels > 0)
        p = zeros(npar, niter);
        p(:,1) = p0(:, 1);
        p_theta = Inf(niter, 1).*-1;
        p_theta(1) = p_theta1;
        for m = 2:niter

            p(:,m) = p(:, find(max(p_theta)));
            p(labels == 1, m) = p(labels == 1, m) + random (pf.noise.type_p, 0, pf.noise.sigma_p(labels == 1))';
            fprintf('Iteration = %d/%d\n',m,niter);
            % initial states estimate
            xh(:, 1) = x0;
            yh(:, 1) = y0;
            for j = 1: nx
                for i = 1: np                          % simulate initial particles
                    xp(j,1,i) = x0(j) + random (pf.noise.type_x, 0, pf.noise.sigma_x(j));
                end
            end
            wxp(1,:) = repmat(1/np, np, 1);
            uwxp(1,:) = repmat(1/np, np, 1); % all particles have the same weight
            wxk = 1/np*ones(np,1);
            for k = 2: T
                %fprintf('Time point = %d/%d\n',k,T);
                % states estimates
                for i = 1: np
                    w = zeros(nx,1);
                    for j = 1: nx
                        w(j) = random (pf.noise.type_x, 0, pf.noise.sigma_x(j));
                    end
                    if isempty(u)
                        xp(:,k,i) = xp(:,k-1,i) + (f(k, p(:,m), xp(:,k-1,i), []) + w) * Dt;
                        yp(:,k,i) = g(k, p(:,m), xp(:,k,i), []);
                    else
                        xp(:,k,i) = xp(:,k-1,i) + (f(k, p(:,m), xp(:,k-1,i), u(:,k)) + w) * Dt;
                        yp(:,k,i) = g(k, p(:,m), xp(:,k,i), u(:,k));
                    end
                    % weights (when using the PRIOR pdf): eq 63, Ref 1
                    wxp(k,i) = wxk(i) * evalpdfobs( y(:,k), yp(:,k,i), pf.noise );  %%%% mean of p's yp given xp

                end
                % Normalize weight vector
                uwxp(k,:) = wxp(k,:);
                wxp(k,:) = wxp(k,:)./sum(wxp(k,:),2);
                wxk = wxp(k,:);
                xk = squeeze(xp(:,k,:));
                Neff = 1/sum(wxk.^2);
                resample_percentage = 0.5;
                Nt  = resample_percentage*np;
                if Neff < Nt
                    % Resampling
                    %disp('States Particles Resampling ...')
                    % remove this condition and sample on each iteration:
                    [xp(:,k,:), wxk,~] = resample(xk, wxk', pf.resampling_strategy);
                end
                % Compute estimated state
                xh(:,k) = squeeze(xp(:,k,:))*wxk';
                if isempty(u)
                    % filtered observation
                    yh(:,k) = g(k, p(:,m), xh(:,k),[]);
                else
                    % filtered observation
                    yh(:,k) = g(k, p(:,m), xh(:,k), u(:,k));
                end
            end
%             p_theta(m) = prod(sum(uwxp,2));
%             p_theta(m) = geomean(sum(uwxp,2));
            p_theta(m) = sum(log((sum(uwxp,2))));
        end

        fprintf('Best State Estimation');
        % initial states estimate
        xh(:, 1) = x0;
        yh(:, 1) = y0;
        for j = 1: nx
            for i = 1: np                          % simulate initial particles
                xp(j,1,i) = x0(j) + random (pf.noise.type_x, 0, pf.noise.sigma_x(j));
            end
        end
        wxp(1,:) = repmat(1/np, np, 1);
        uwxp(1,:) = repmat(1/np, np, 1); % all particles have the same weight
        wxk=1/np*ones(np,1);

        for k=2:T
            fprintf('Time point = %d/%d\n',k,T);
            % states estimates
            for i = 1:np
                w=zeros(nx,1);
                for j=1:nx
                    w(j) = random (pf.noise.type_x, 0, pf.noise.sigma_x(j));
                end
                if isempty(u)
                    xp(:,k,i) = xp(:,k-1,i) + (f(k,p(:, find(max(p_theta))), xp(:,k-1,i), []) + w) * Dt;
                    yp(:,k,i) = g(k,p(:, find(max(p_theta))), xp(:,k,i), []);
                else
                    xp(:,k,i) = xp(:,k-1,i) + (f(k,p(:, find(max(p_theta))), xp(:,k-1,i), u(:,k)) + w) * Dt;
                    yp(:,k,i) = g(k,p(:, find(max(p_theta))), xp(:,k,i), u(:,k));
                end
                % weights (when using the PRIOR pdf): eq 63, Ref 1
                wxp(k,i) = wxk(i) * evalpdfobs( y(:,k), yp(:,k,i), pf.noise );  %%%% mean of p's yp given x
            end
            % Normalize weight vector
            uwxp(k,:) = wxp(k,:);
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
                yh(:,k) = g(k, p(:, find(max(p_theta))), xh(:,k),[]);
            else
                % filtered observation
                yh(:,k) = g(k, p(:, find(max(p_theta))), xh(:,k),u(:,k));
            end
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