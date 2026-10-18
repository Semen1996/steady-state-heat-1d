function [T, iter, residual] = gauss_seidel(T, tol, maxIterations)
  % [T, iter, residual] = gauss_seidel(T, tol, maxIterations)
  %
  % Gauss-Seidel method for the system T(i-1) - 2T(i) + T(i+1) = 0,
  % i = 2..nx-1. Nodes are swept from left to right, and the already
  % updated value of the left neighbour is used immediately.
  %
  % Input:
  %   T             - row vector of length nx: initial guess;
  %                   T(1) and T(end) are the boundary temperatures TL and TR
  %   tol           - threshold for max|T^(k+1) - T^(k)| (default 1e-5)
  %   maxIterations - maximum number of iterations (default 10000)
  %
  % Output:
  %   T        - solution (boundary values are unchanged)
  %   iter     - number of iterations performed
  %   residual - residual max|T(i-1) - 2T(i) + T(i+1)| over interior nodes
  %
  % Issues a warning if the criterion max|T^(k+1) - T^(k)| <= tol
  % is not met within maxIterations iterations.


  if(nargin < 2 || isempty(tol))
    tol = 1e-5;
  end
  if(nargin < 3 || isempty(maxIterations))
    maxIterations = 10000;
  end
  nx = length(T);

  TOld = T;
  converged = false;

  for iter = 1:maxIterations
    for i = 2:nx-1;
      T(i) = (T(i+1) + T(i-1)) /2;
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
    warning('gauss_seidel: iteration limit reached (%d), max|dT| = %g, residual = %g', ...
            maxIterations, delta, residual);
  end
end

