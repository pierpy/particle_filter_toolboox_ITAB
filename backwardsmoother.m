function [ xhs, yhs, xps, wps, uwps ] = backwardsmoother( t, p0, u, f, wxp, xp, xh, yh, pf )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    T = length(xh);
    Dt = mean(diff(t));  
    np = size(wxp,2);
    nx = size(xp,1);
    xhs = zeros(size(xh));
    xps = zeros(size(xp));
    wps = zeros(nx, size(wxp,1), size(wxp,2));
    yhs = zeros(size(yh));
    xhs(:, T) = xh(:, T);
    for i = 1: nx
         wps(i, T, :) = wxp(T, :);
    end
    for tt = T-1:-1:1    
%         fprintf('Time point = %d/%d\n',tt,T);
        for jj = 1: np
            sev = squeeze(xp(:, tt, :) + (f(tt, p0, xp(:, tt, :), u(tt))) * Dt);
            for kk = 1:nx
                tmp = normpdf(squeeze(xp(kk, tt+1, jj)), sev(kk, :)', pf.noise.sigma_x);
                PDF(kk,:) = tmp(:,1)';
            end
            v(jj,tt) = sum(wxp(tt,:).*mean(PDF));
        end
        
        % Compute the smoothing weight
        for ii = 1: np
            sev = xp(:, tt, ii) + (f(tt, p0, xp(:, tt, ii), u(tt))) * Dt;
            for kk = 1: nx
                tmp = normpdf(squeeze(xp(kk, tt+1, :)), sev(kk), pf.noise.sigma_x);
                PDF2(kk, :) = tmp(:, 1)'; 
                wps(kk, tt, ii) = wxp(tt,ii)*sum(squeeze(wps(kk, tt+1, :))'.*PDF2(kk, :)./v(:,tt)');
            end       
%             wps(tt,ii) = wxp(tt,ii)*sum(wps(tt+1, :).*mean(PDF)./v(:,tt)');
        end
        % normalization term
%         for i = 1: np
%             v(jj,tt)=sum(wxp(tt,:).*normpdf(xp(tt+1, i),sys.a*p(:,tt),sys.sigmav));  
%         end
%         for j = 1: np
%             sev = xp(:,tt,j) + (f(tt, p0, xp(:, tt, j), u(tt))) * Dt;
%             wps(tt, j) = wxp(tt, j)*mean(normpdf(xp(tt+1,:), sev, pf.noise.sigma_x));
%         end
        uwps = wps;
        for i =1: nx
            wps(i, tt, :) = wps(i, tt, :)./sum(wps(i, tt, :));       
            xhs(i, tt) = squeeze(wps(i, tt, :))' * squeeze(xp(i, tt, :));
        end        
    end
end

