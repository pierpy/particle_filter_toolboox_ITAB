% Particle filter

function str = ppff(m, M, z, u, f, g, theta)
  y           = z.y;
  N           = size(y,2);   % Number of samples
  nx = size(m.X1,1);
  ny = size(y, 1);
  xf          = zeros(nx, size(y,2));
  xPWeighted  = zeros(nx, N, M);      % Make some room
  xPResampled = zeros(nx, N, M);      % Make some room
  wStore      = zeros(N, M);      % Make some room

  % 1. Initialize particles
  if m.P1==0
    x = repmat(m.X1,1,M);    % Initial state fully known
  else
    x = repmat(m.X1,1,M) + m.Q(1)*randn(nx, M);
  end;
  for t = 1: N
    for i = 1: M
        e(:, i) = y(:, t) - g( t, theta, x(:, i), u(:, t) ); 
        w(i) = exp(-(1/2)*(e(:, i)'*inv(m.R)*e(:, i)));
    end% 2. Evaluate the likelihood
%     w = mean(exp(-(1/2)*(e.*(m.R\e))), 1);  % 2. Evaluate the likelihood
    w = w/sum(w);
    xPWeighted(:, t, :) = x;
    wStore(t,:)  = w;
    xf(:, t) = (w*x')';             % Compute state estimate
    index = sysresample1(w);        % 3. Resample
    for i = 1: size(x, 1)
        x(i, :) = x(i, index);
    end
    xPResampled(:, t, :) = x;

        for i = 1: M
    %         x(:, i) = x(:, i) + f( t, theta, x(:, i), u(:, t) ) + sqrt(m.Q)*randn(m.nstates, 1); % 4. Predict the particles
              x(:, i) = x(:, i) + (f( t, theta, x(:, i), u(:, t) ) + m.Q(1).*randn(nx,1))*m.dt;
%             x(:, i) = x(:, i) + (f( t, theta, x(:, i), u(:, t) ) +sqrt(m.Q)*randn(nx,1))*m.dt;
        end

  end;
  str.Xf          = xf;
  str.xPWeighted  = xPWeighted;
  str.xPResampled = xPResampled;
  str.w           = wStore;
end
