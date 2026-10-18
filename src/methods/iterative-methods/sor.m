function [T, iter, residual] = sor(T, tol, maxIterations, omega)
  % [T, iter, residual] = sor(T, tol, maxIterations, omega)
  %
  % Successive over-relaxation (SOR) method for the system
  % T(i-1) - 2T(i) + T(i+1) = 0, i = 2..nx-1.
  % The Gauss-Seidel correction is multiplied by omega; omega = 1
  % gives the Gauss-Seidel method.
  %
  % Input:
  %   T             - row vector of length nx: initial guess;
  %                   T(1) and T(end) are the boundary temperatures TL and TR
  %   tol           - threshold for max|T^(k+1) - T^(k)| (default 1e-5)
  %   maxIterations - maximum number of iterations (default 10000)
  %   omega         - relaxation parameter, 0 < omega < 2 (default 1.5)
  %
  % Output:
  %   T        - solution (boundary values are unchanged)
  %   iter     - number of iterations performed
  %   residual - residual max|T(i-1) - 2T(i) + T(i+1)| over interior nodes
  %
  % Issues a warning if the criterion max|T^(k+1) - T^(k)| <= tol
  % is not met within maxIterations iterations.


  nx = length(T);
  if(nargin < 2 || isempty(tol))
    tol = 1e-5;
  end
  if(nargin < 3 || isempty(maxIterations))
    maxIterations = 10000;
  end
  if(nargin < 4 || isempty(omega))
    omega = 1.5;
  end

  TOld = T;
  converged = false;

  for iter = 1:maxIterations
    for i = 2:nx-1;
      TStar = (T(i+1) + T(i-1)) / 2;
      T(i) = T(i) + omega*(TStar - T(i));
    end

    delta = max(abs(TOld - T), [], 'all');

    if(delta <= tol)
      converged = true;
      break;
    end

    TOld = T;
  end

  residual = max(abs(T(1:end-2) - 2*T(2:end-1) + T(3:end)));

  if(~converged)
    warning('sor: iteration limit reached (%d), max|dT| = %g, residual = %g', ...
            maxIterations, delta, residual);
  end
end

