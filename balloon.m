function f=balloon(t, p, x, u)
    % initial states
    x1 = x(1);                               % vasodilatory signal s(t)
    x2 = x(2);                               % blood inflow f(t)
    x3 = x(3);                               % blood volume v(t)
    x4 = x(4);
    x5 =  x(5);
    fv = x3.^(1./p(5,1));                    % blood outflow                   
    ff = (1-(1-p(6,1)).^(1./x2))./p(6,1);    % oxygen extraction fraction   
    f = zeros(5,1);
    % Evaluate flow field
    f(1) = (p(1,1).*x5 - p(2,1).*x1 - p(3,1).*((x2) - 1));    % vasodilatatory signal
    f(2) = x1./(x2);                                          % local normalized blood flow
    f(3) =(x2 - fv)./(p(4,1).*(x3));                          % HbO normalized
    f(4) =(x2.*ff./(x4) - fv./(x3))./p(4,1);                  % Hbb normalized
    f(5) = p(7,1)*x5;                                         % neuronal activity
end