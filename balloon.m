function f=balloon(t,p,x,u)


%% parameters
lambda=-0.3; % autocorrelation neuronal
epsilon=1; % neuronal efficacy
kas=0.2632; % signal decay
kaf=0.1297; %autoregulation
tau0=3.14; %transit time
alpha=0.3663; %stiffness
E0=0.4865; % resting Oxygen extraction

    x1 = x(1);                               % vasodilatory signal s(t)
    x2 = x(2);                     % blood inflow f(t)
    x3 = x(3);                     % blood volume v(t)
    x4 = x(4);
    x5 =  x(5);
   
fv =x3.^(1./alpha);                            % blood outflow                   
ff = (1-(1-E0).^(1./x2))./E0;                % oxygen extraction fraction   
f = zeros(5,1);
% Evaluate flow field
f(1) = (epsilon.*x5 - kas.*x1 - kaf.*((x2) - 1));         % vasodilatatory signal
f(2) = x1./(x2);                                          % local normalized blood flow
f(3) =(x2 - fv)./(tau0.*(x3));                              % HbO normalized
f(4) =(x2.*ff./(x4) - fv./(x3))./tau0;                        % Hbb normalized
f(5)=lambda*x5+randn/10^2;                                    % neuronal activity


end