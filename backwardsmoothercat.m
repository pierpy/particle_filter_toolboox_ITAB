function [ xps, wps, xhs ] = backwardsmoothercat( t, p0, u, f, g, wxp, xp, xh, yh, y, pf )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    T = length(xh);
    Dt = mean(diff(t));  
    np = size(wxp,2);
    nx = size(xp,1);
    xhs = zeros(size(xh));
    xps = zeros(size(xp));
    wps = zeros(size(wxp));
    yhs = zeros(size(yh));
    xhs(:, T) = xh(:, T);
    wps(T, :) = wxp(T, :);

    multinomial1 = makedist('Multinomial', 'probabilities', wxp(T, :));
    wind = random(multinomial1, 1, np);
    for i = 1: nx
        xps(i, T, :) = xp(i, T, wind);
        
    end
    for tt = T-1:-1:1    
%         fprintf('Time point = %d/%d\n',tt,T);
        for jj = 1: np
            multinomial = makedist('Multinomial', 'probabilities', wxp(tt, :));
            a = random(multinomial);  
            for ii = 1: np
                updatedWeights(ii) = mean(wps(tt+1, ii).* g(tt+1, p0, xps(:, tt+1, ii), u(tt)));
            end     
            updatedWeights = updatedWeights./sum(updatedWeights);
            multinomialw = makedist('Multinomial', 'probabilities', updatedWeights);
            b = random(multinomialw); 
            xps(:, tt, jj) = xp(:, tt, a); 
            sev = xps(:, tt+1, b) + (f(tt, p0, xps(:, tt, jj), u(tt))) * Dt;
            for i = 1:nx
                pp(i) = normpdf(sev(i), xp(i, tt, jj), pf. noise.sigma_x(i));
            end
            wps(tt, jj) = mean(pp);
        end
      wps(tt, :) = wps(tt, :)./sum(wps(tt, :));
    end
    for t = 1:T
        xhs(:, t) = wps(t, :)*squeeze(xps(:, t, :));
    end
end

