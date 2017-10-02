function [ xh, yh, uwxp, wxp, xp, xpnr ] = particlefilter( t, x0, p0, u, y, f, g, pf )
%UNTITLED3 Summary of this function goes here
    %   Detailed explanation goes here
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
%         fprintf('Time point = %d/%d\n',k,T);
        % states estimates
        for i = 1: np
            w = zeros(nx,1);
            for j=1:nx
                w(j) = random (pf.noise.type_x, 0, pf.noise.sigma_x(j));
            end
            if isempty(u)
                xp(:, k, i) = xp(:, k-1, i) + (f(k, p0, xp(:,k-1,i), []) + w) * Dt;
                yp(:, k, i) = g(k, p0, xp(:, k, i), []);
            else
                xp(:, k, i) = xp(:, k-1, i) + (f(k, p0, xp(:,k-1,i), u(:,k)) + w) * Dt;
                yp(:, k, i) = g(k, p0, xp(:, k, i), u(:,k));
            end
            % weights (when using the PRIOR pdf): eq 63, Ref 1
            wxp(k,i) = wxk(i) * evalpdfobs( y(:, k), yp(:, k, i), pf.noise );  %%%% mean of p's yp given xp
            
        end
        % Normalize weight vector
        uwxp(k,:) = wxp(k,:);
        wxp(k,:) = wxp(k,:)./sum(wxp(k,:),2);
        wxk = wxp(k,:);
        xpnr = xp;
        xk = squeeze(xp(:,k,:));        
        Neff = 1/sum(wxk.^2);
        resample_percentage = 0.5;
        Nt  = resample_percentage*np;
        if Neff < Nt
            % Resampling
%             disp('States Particles Resampling ...')
            % remove this condition and sample on each iteration:
            [xp(:,k,:), wxk, ~] = resample(xk, wxk', pf.resampling_strategy);
        end
        % Compute estimated state
        if nx == 1
          xh(:,k) = squeeze(xp(:,k,:))'*wxk';  
        else
          xh(:,k) = squeeze(xp(:,k,:))*wxk';
        end
        if isempty(u)
            % filtered observation
            yh(:,k) = g(k, p0, xh(:,k), []);
        else
            % filtered observation
            yh(:,k) = g(k, p0, xh(:,k), u(:,k));
        end
    end
%     loglikelihood = sum(log((sum(uwxp,2))));
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
if size(xk,2)==1
    xk = xk(idx,:);
    wk = repmat(1/Ns,1,Ns);
else
    xk = xk(:,idx);
    wk = repmat(1/Ns,1,Ns);
end

end
