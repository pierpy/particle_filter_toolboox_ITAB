function f=balloonInput(t,p,x,u)


%% parameters


    x1 = x(1);                               % vasodilatory signal s(t)
    x2 = x(2);                     % blood inflow f(t)
    x3 = x(3);                     % blood volume v(t)
    x4 = x(4);
    
   
fv = x3.^(1./p(5));                            % blood outflow                   
ff = (1-(1-p(6)).^(1./x2))./p(6);                % oxygen extraction fraction   
f = zeros(4,1);
% Evaluate flow field
    f(1) = (p(1).*u - p(2).*x1 - p(3).*((x2) - 1));    % vasodilatatory signal
    f(2) = x1./(x2);                                          % local normalized blood flow
    f(3) =(x2 - fv)./(p(4).*(x3));                          % HbO normalized
    f(4) =(x2.*ff./(x4) - fv./(x3))./p(4);   % Hbb normalized
%     f(4) = x2./(p(4,1)*x4)*(ff - fv./x3);
                                 


end