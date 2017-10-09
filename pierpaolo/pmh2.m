function [ loglikelihood, theta, proposedThetaAccepted ] = pmh2( initialTheta, model)
% PMH Summary of this function goes here
    %   Detailed explanation goes here
        theta = zeros(size(initialTheta, 1), model.niter); 
        thetaProposed = zeros(size(initialTheta, 1), model.niter);
        proposedThetaAccepted = zeros(size(initialTheta, 1), model.niter);
        loglikelihood = zeros(1, model.niter);
        loglikelihoodProposed = zeros(1, model.niter);
        theta(:, 1) = initialTheta;

        [ model ] = pFilter( model.t, model );
        loglikelihood(1) = model.filter.filtLoglikelihood;
        for m = 2: model.niter % run MPH
            tic
            thetaProposed(:, m) = theta(:, m-1);
            thetaProposed(model.labels == 1, m) = mvnrnd(theta(model.labels == 1, m-1), model.stepSize); 
            
            fprintf('Iteration = %d/%d\n', m, model.niter);
            
            if abs(thetaProposed(model.labels == 1, m)) > 2.9
                [ model ] = pFilter( model.t, model );
                loglikelihoodProposed(m) = model.filter.filtLoglikelihood;
            end
            
            prior = dnorm(thetaProposed(model.labels == 1, m), 3, 0.01);
            prior = prior - dnorm(theta(model.labels == 1, m-1), 3, 0.01);
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