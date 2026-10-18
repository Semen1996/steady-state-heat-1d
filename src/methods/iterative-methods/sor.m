function [U, iter, residual] = sor(U, tol, maxIterations, omega)
  % iter     - число выполненных итераций
  % residual - max|T(i-1) - 2T(i) + T(i+1)| по внутренним узлам
  % tol           - порог max|T^(k+1) - T^(k)| (по умолчанию 1e-5)
  % maxIterations - наибольшее число итераций (по умолчанию 10000)
  % omega         - параметр релаксации, 0 < omega < 2 (по умолчанию 1.5)
  nx = length(U);
  if(nargin < 2 || isempty(tol))
    tol = 1e-5;
  end
  if(nargin < 3 || isempty(maxIterations))
    maxIterations = 10000;
  end
  if(nargin < 4 || isempty(omega))
    omega = 1.5;
  end

  UOld = U;
  converged = false;

  for iter = 1:maxIterations
    for i = 2:nx-1;
      UStar = (U(i+1) + U(i-1)) / 2;
      U(i) = U(i) + omega*(UStar - U(i));
    end

    delta = max(abs(UOld - U), [], 'all');

    if(delta <= tol)
      converged = true;
      break;
    end

    UOld = U;
  end

  residual = max(abs(U(1:end-2) - 2*U(2:end-1) + U(3:end)));

  if(~converged)
    warning('sor: достигнут лимит итераций (%d), max|dT| = %g, невязка = %g', ...
            maxIterations, delta, residual);
  end
end

