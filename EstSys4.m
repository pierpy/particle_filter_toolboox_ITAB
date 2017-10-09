%% Particle Filter State Space Model Evaluation with 
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
% niter=number of parameter estimation
% labels = length of p0 vector with 1 estimate 0 do not estimate


function [xh, yh, ph, loglikelihood] = EstSys4 (t, x0, p0, p, labels, niter, u, y, f, g, pf)

    
    
    % initial states estimate
    [ ~, ~, uwxp, loglikelihood1 ] = particlefilter(t, x0, p0, u, y, f, g, pf );
    if sum(labels > 0)
    [ loglikelihood, ph] = pmh( t, x0, p0, p, labels, niter, u, y, f, g, pf, uwxp, loglikelihood1);
    
    [ xh, yh, ~, ~ ] = particlefilter(t, x0, ph(:, find(max(loglikelihood))), u, y, f, g, pf );

    end
end
