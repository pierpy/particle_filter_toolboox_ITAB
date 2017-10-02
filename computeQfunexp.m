% Used in solving the nonlinear identification problem considered in 
% Section 10.2 in
%
% Thomas B. Schön, Adrian Wills and Brett Ninness. System Identification 
% of Nonlinear State-Space Models. Automatica, 47(1):39-49, January 2011. 
%
% theta  = current value for the parameters.
% xPart  = particles to be used in order to form the Q function.
% y      = measurements to be used in order to form the Q function.
% 
% Qfun   = the value of the Q function at the point theta
%


% x(t+1) = m.a*x(t) + m.b*x(t)/(1+x(t)^2) + m.c*cos(1.2*t) + v(t);
%   y(t)   = m.d*x(t)^2 + e(t);

function Qfun = computeQfunexp(p0, xps, y, f, g, u, pf)

  Qx = diag(pf.noise.sigma_x);
  Qy = diag(pf.noise.sigma_y);
  N = size(xps, 2);   % Number of samples
  M = size(xps, 3);   % Number of particles 
  np = size(xps, 3);
  nx = size(xps, 1);
  Qfun1 = 0;
  Qfun2 = 0;      
  for t = 1: N-1
      for p = 1: np
          if isempty (u)
              Qfun1tmp(:, p) = (squeeze(xps(:, t+1, p)) - f(t, p0, squeeze(xps(:, t, p)),[]) - xps(:, t, p)).^2;
          else
        Qfun1tmp(:, p) = (squeeze(xps(:, t+1, p)) - f(t, p0, squeeze(xps(:, t, p)), u(:,t)) - xps(:, t, p)).^2;     
          end
      end
      Qfun1tmp = sum(Qfun1tmp);
      Qfun1    = Qfun1 + (1/(2*M*Qx*Qx))*Qfun1tmp;
  end
  for t = 1: N
      for p = 1: np
          if isempty(u)
          Qfun2tmp = (y(t) - g(t, p0, xps(:, t, p), [])).^2;
          else
          Qfun2tmp = (y(t) - g(t, p0, xps(:, t, p), u(:,t))).^2;
          end
      end
    Qfun2tmp = sum(Qfun2tmp);
    Qfun2    = Qfun2 + (1/(2*M*Qy*Qy))*Qfun2tmp;
  end;
  Qfun = -((N-1)/2)*log(Qx*Qx) - Qfun1 - (N/2)*log(Qy*Qy) - Qfun2;
  Qfun = -Qfun;       % We want to maximize Q, but fminunc solves a minimization problem
end