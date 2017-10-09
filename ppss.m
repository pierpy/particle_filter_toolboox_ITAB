% Particle smoother based on Douc, Garivier, Moulines and Olsson (2010)

function str = ppss(Q, M, particles, u, f, g, theta)
  N     = size(particles,2);            % Number of samples
  Qi    = inv(Q);
%  B     = 1/sqrt(((2*pi)^m.nx)*det(m.ss.Q));  % Upper bound for the pdf of the dynamics
%  const = 1/sqrt(((2*pi)^m.nx)*det(m.ss.Q));
  J     = zeros(N, M);              % The index used for the smoothed estimate
  xs    = zeros(size(particles,1), N);
  xPart = zeros(size(particles,1), N, M);
  nrR   = zeros(N, M);              % Keep track of how many rejections we get

  J(N,:) = floor(rand(1, M)*M) + 1;   % Randomly select one of the equaly likely indecies
  for t = N-1:-1:1

    for i = 1: M

      notfound   = 1;
      nrRejected = 0;
      while notfound
        U      = rand(1);
        I      = ceil(rand(1)*M);     % Randomly select one of the equaly likely indecies
%        w      = particles(t+1,J(t+1,i)) - m.ss.A*particles(t,I);
        xt     = particles(:, t, I);
        xt1    = particles(:, t+1, J(t+1, i));
        w      = xt1 - xt - f( t, theta, xt, u(:, t) ) ;
        target = exp(-0.5*w'*Qi*w);   % Evaluate the target density for the choosen index
        if U <= target
          J(t,i)   = I;
          W(t, i) =target;
          notfound = 0;
        else
          nrRejected = nrRejected + 1;
        end;
      end;
      nrR(t,i) = nrRejected;
    end;
  end;

  for t = 1: N
    xPart(:, t, :) = particles(:, t,  J(t,:));
    xPart1 = squeeze(particles(:, t,  J(t,:)));
    if size(xs, 1) == 1
        xs(:, t)      = mean(xPart1);
    else
        xs(:, t)      = mean(xPart1, 2);
    end
  end;
  str.xPart = xPart;      % Return the particles
  str.Xs    = xs;         % Return the smoothed states estimate
  str.J     = J;          % Return the resulting index
  str.nrR   = nrR; 
  str.w = W;% Return the number of rejections
end
