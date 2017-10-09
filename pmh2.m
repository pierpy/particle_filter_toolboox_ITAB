function [ loglikelihood, theta, proposedThetaAccepted ] = pmh2( t, x0, initialTheta, labels, niter, u, y, f, g, pf, stepSize )
% PMH Summary of this function goes here
    %   Detailed explanation goes here
        theta = zeros(size(initialTheta, 1), niter); 
        thetaProposed = zeros(size(initialTheta, 1), niter);
        proposedThetaAccepted = zeros(size(initialTheta, 1), niter);
        loglikelihood = zeros(1, niter);
        loglikelihoodProposed = zeros(1, niter);
        theta(:, 1) = initialTheta;

        [ ~, loglikelihood(1) ] = particlefilter( t, x0, theta(:, 1), u, y, f, g, pf );
        for m = 2: niter % run MPH
            tic
            thetaProposed(:, m) = theta(:, m-1);
            thetaProposed(labels == 1, m) = mvnrnd(theta(labels == 1, m-1), stepSize); 
            
            fprintf('Iteration = %d/%d\n', m, niter);
            
            if abs(thetaProposed(labels == 1, m)) > 2.9
                [ ~ , loglikelihoodProposed(m) ]  = particlefilter( t, x0, thetaProposed(:, m), u, y, f, g, pf );
            end
            
            prior = dnorm(thetaProposed(labels == 1, m), 3.3, 0.1);
            prior = prior - dnorm(theta(labels == 1, m-1), 3.3, 0.1);
            likelihoodDifference = loglikelihoodProposed(m) - loglikelihood(m - 1);
            acceptProbability = exp(prior + likelihoodDifference); 
            acceptProbability = acceptProbability * (abs(thetaProposed(4, m)) < 4); 
            uniformRandomVariable = unifrnd(0, 1);
            
            if (uniformRandomVariable < acceptProbability)
                % Accept the parameter
                theta(:, m) = thetaProposed(:, m);
                loglikelihood(m) = loglikelihoodProposed(m);
                proposedThetaAccepted(m) = 1.0;
            else
                % Reject the parameter
                theta(:, m) = theta(:, m-1);
                loglikelihood(m) = loglikelihood(m - 1);
                proposedThetaAccepted(m) = 0.0;
            end  
            toc
        end       
end

% helper functions
function[out] = dnorm(x, mu, sigma)
    out = -0.5 .* log(2 * pi) - 0.5 .* log(sigma.^2) - 0.5 ./ sigma.^2 .* (x - mu).^2;
end