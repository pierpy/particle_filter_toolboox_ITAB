function [ loglikelihood, ph ] = pmh( t, x0, p0, labels, niter, u, y, f, g, pf )
% PMH Summary of this function goes here
    %   Detailed explanation goes here
        ph = zeros(size(p0,1), niter);
        ph(:,1) = p0;
        loglikelihood = -Inf(niter, 1);
        [ ~, ~, uwxp] = particlefilter(t, x0, p0, u, y, f, g, pf );
        loglikelihood(1) = sum(log((sum(uwxp,2))));
        for m = 2: niter
            ph(:, m) = ph(:, m-1);
            ph(labels == 1, m) = ph(labels == 1, m-1) + random (pf.noise.type_p, 0, pf.noise.sigma_p(labels == 1))';
            fprintf('Iteration = %d/%d\n', m, niter);
            % initial states estimate
            [ ~, ~, uwxp] = particlefilter( t, x0, ph(:, m), u, y, f, g, pf );
            loglikelihood(m) = sum(log((sum(uwxp,2))));
            al = min(0, (loglikelihood(m) - loglikelihood(m-1)));
            uu = log(rand);
            if uu < al && isreal(loglikelihood(m))
                ph(:, m) = ph(:, m);   
                disp(strcat('uu < al ', num2str(loglikelihood(m)), ' \', num2str(ph(:,m))))
            else
                ph(:, m) = ph(:, m-1);
                loglikelihood(m) = loglikelihood(m-1);
            end
        end       
end

