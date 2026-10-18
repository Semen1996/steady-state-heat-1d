function [U, iter, residual] = sor(U)
  % iter     - число выполненных итераций
  % residual - max|T(i-1) - 2T(i) + T(i+1)| по внутренним узлам
  nx = length(U);
  maxIterations = 10000;
  omega = 1.5;

  UOld = U;
  converged = false;

  for iter = 1:maxIterations
    for i = 2:nx-1;
      UStar = (U(i+1) + U(i-1)) / 2;
      U(i) = U(i) + omega*(UStar - U(i));
    end

    delta = max(abs(UOld - U), [], 'all');

    if(delta <= 1e-5)
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

