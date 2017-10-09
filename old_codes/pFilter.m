function [ xh, loglikelihood, particles ] = pFilter( t, x0, p0, u, y, f, g, pf )
%UNTITLED3 Summary of this function goes here
    %   Detailed explanation goes here
    model.T = length(t);
                              % signal length
    model.Dt = mean(diff(t));               % sampling time
    model.nx = length(x0);                  % number of states
    if isempty(u)
        y0 = g(t(1), p0(:,1), x0, []);           % initial observation
    else
        y0 = g(t(1), p0(:,1), x0, u(:,1));    % initial observation
    end
    model.ny = length(y0);                  % number of observation
    model.npar = length(p0(:,1));                % number of parameters
    model.np = pf.np; % number of particles
    model.p0 = p0;
    model.f = f;
    model.g = g;
    model.pf = pf;
    model.Q = pf.Q;
    model.R = pf.R;
    model.u = u;
    % initialization of state, observation and pf.noise
    xh = zeros(model.nx, model.T);
    yh = zeros(model.ny, model.T);
    ph = zeros(model.npar, model.T);
    xp = zeros(model.nx, model.np);
    yp = zeros(model.ny, model.np);
    wxp = zeros(model.T, model.np);
    uwxp = zeros(model.T, model.np);
    
    particles = zeros(model.nx, model.T, model.np);

    xh(:, 1) = x0;
    yh(:, 1) = y0;
    loglikelihood = 0;
                        % simulate initial particles
    xp = repmat(x0', 1, model.np) + pf.Q*randn(model.nx, model.np);
    particles(:, 1, :) = xp;

    wxp(1,:) = repmat(1/model.np, model.np, 1); 
    uwxp(1,:) = repmat(1/model.np, model.np, 1); % all particles have the same weight
    wxk = 1/model.np * ones(model.np, 1);
    for k=2: model.T
%         fprintf('Time point = %d/%d\n',k,T);
        % states estimates     
        xp = predictParticles(k, xp, model);
        yp = measUpdate(k, xp, model);
        wxp(k, :) = computeWeights(y(:, k), yp, model.R);
        % Normalize weight vector
        maxw = max(wxp(k, :));
        wxp(k,:) = wxp(k, :)./sum(wxp(k, :),2);
        wxk = wxp(k,:);
        xk = xp;        
        Neff = 1/sum(wxk.^2);
        resample_percentage = 0.5;
        Nt  = resample_percentage*model.np;
        if Neff < Nt
            % Resampling
%             disp('States Particles Resampling ...')
            % remove this condition and sample on each iteration:
            [xp, wxk, ~] = resample(xk, wxk', pf.resampling_strategy);
        end
        particles(:, k, :) = xp;
        % Compute estimated state
        if model.nx == 1
          xh(:,k) = xp*wxk';  
        else
          xh(:,k) = xp*wxk';
        end
        if isempty(model.u)
            % filtered observation
            yh(:,k) = g(k, p0, xh(:,k), []);
        else
            % filtered observation
            yh(:,k) = g(k, p0, xh(:,k), u(:,k));
        end
        predectivelikelihood = maxw + log(sum(wxp(k,:),2)) - log(model.np);
        loglikelihood = loglikelihood + predectivelikelihood;
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
