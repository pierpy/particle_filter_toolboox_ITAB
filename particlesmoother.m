function [ J, xps, xhs] = particlesmoother( t, p0, u, f, xp, pf )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    T = length(t);
    np = size(xp, 3);
    nx = size(xp, 1);
    % creo la matrice delle covariaze che mi servirà per calcolo pdf
    Qx = diag(pf.noise.sigma_x);
    Qy = diag(pf.noise.sigma_y);
    Qp = diag(pf.noise.sigma_p);    
    J     = zeros(T,np);              % The index used for the smoothed estimate
    nrR   = zeros(T,np);              % Keep track of how many rejections we get
    J(T,:) = floor(rand(1,np)*np)+1;  % initial indexes
    
    xhs = zeros(size(xp,1),T);
    xps = zeros(size(xp));

    for t = T-1:-1:1
%         fprintf('Time point = %d/%d\n',t,T);
        for p = 1: np
            disp(['particle:', num2str(p)])
            notfound   = 1;
            nrRejected = 0;
            while notfound
                disp(['notfound:', num2str(p)])
                U = rand(1);
                I = ceil(rand(1)*np);
                xt = xp(:, t, I);
                xt1 = xp(:, t+1, J(t+1, p));
                if isempty(u)
                    w = xt1 - f(t, p0, xt, []) - xt;
                else
                w = xt1 - f(t, p0, xt, u(:,t)) - xt;
                end
%                 target = mvnpdf(w', 0, Qx);
                target = exp(-0.5*w'*Qx*w); 
                if U <= target
                    J(t, p)   = I;
                    notfound = 0;
                else
                    nrRejected = nrRejected + 1;
                end
            end
            nrR(t, p) = nrRejected;
        end
    end
    for s = 1: nx
        for t = 1: T
            xps(s, t, :) = xp(s, t, J(t,:));
            xhs(s, t) = mean(xps(s, t,:));
        end
    end
end

