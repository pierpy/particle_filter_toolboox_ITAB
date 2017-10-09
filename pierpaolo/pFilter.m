function [ model ] = pFilter( t, model )
%UNTITLED3 Summary of this function goes here
    %   Detailed explanation goes here
                     % number of states
    if isempty(model.u)
        y0 = model.g(t(1), model.p0(:,1), model.x0, []);           % initial observation
    else
        y0 = model.g(t(1), model.p0(:,1), model.x0, model.u(:,1));    % initial observation
    end
    model.ny = length(y0);  
    model.y0 = y0;
    % initialization of state, observation and pf.noise
    xh = zeros(model.nx, model.T);
    yh = zeros(model.ny, model.T);
    ph = zeros(model.npar, model.T);
    xp = zeros(model.nx, model.np);
    yp = zeros(model.ny, model.np);
    wxp = zeros(model.T, model.np);
    uwxp = zeros(model.T, model.np);
    
    particles = zeros(model.nx, model.T, model.np);

    xh(:, 1) = model.x0;
    yh(:, 1) = model.y0;
    loglikelihood = 0;
                        % simulate initial particles
    xp = repmat(model.x0', 1, model.np) + model.Q*randn(model.nx, model.np);
    particles(:, 1, :) = xp;

    wxp(1,:) = repmat(1/model.np, model.np, 1); 
    uwxp(1,:) = repmat(1/model.np, model.np, 1); % all particles have the same weight
    wxk = 1/model.np * ones(model.np, 1);
    for k=2: model.T
%         fprintf('Time point = %d/%d\n',k,T);
        % states estimates     
        xp = predictParticles(k, xp, model);
        yp = measUpdate(k, xp, model);
        wxp(k, :) = computeWeights(model.y(:, k), yp, model.R);
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
            [xp, wxk, ~] = resample(xk, wxk', model.resampling_strategy);
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
            yh(:,k) = model.g(k, model.p0, xh(:,k), []);
        else
            % filtered observation
            yh(:,k) = model.g(k, model.p0, xh(:,k), model.u(:,k));
        end
        predectivelikelihood = maxw + log(sum(wxp(k,:),2)) - log(model.np);
        loglikelihood = loglikelihood + predectivelikelihood;
    end
    model.filter.xFilt = xh;
    model.filter.filtLoglikelihood = loglikelihood;
    model.filter.filtParticles = particles;
end
% helper functions
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
function [ predictedParticles ] = predictParticles(k, pastParticles, model)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    w = model.Q*randn(size(pastParticles,1), size(pastParticles,2));
    for i = 1: size(pastParticles, 2)
        if isempty(model.u)
            predictedParticles(:, i) = pastParticles(:, i) + (model.f(k, model.p0, pastParticles(:, i), []) + w(:, i)) * model.Dt;
        else
            predictedParticles(:, i) = pastParticles(:, i) + (model.f(k, model.p0, pastParticles(:, i), model.u(:,k)) + w(:, i)) * model.Dt;
        end
    end

end
function [ yp ] = measUpdate(k, predictedParticles, model )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    for i = 1: size(predictedParticles, 2)
        if isempty(model.u)
            yp(:, i) = model.g(k, model.p0, predictedParticles(:, i), []);
        else
            yp(:, i) = model.g(k, model.p0, predictedParticles(:, i), model.u(:,k));
        end
    end
end
function [ ww ] = computeWeights( y, yp, R )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    for i = 1: size(yp, 2)
        ww(i) = mvnpdf(y - yp(:,i), 0,  R) + 1e-300;
    end

end
