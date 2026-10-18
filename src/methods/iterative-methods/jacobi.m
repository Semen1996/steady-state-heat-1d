function [U, iter, residual] = jacobi(U)
  % iter     - число выполненных итераций
  % residual - max|T(i-1) - 2T(i) + T(i+1)| по внутренним узлам
  nx = length(U);
  maxIterations = 10000;


  UNew = U;
  converged = false;

  for iter = 1:maxIterations
    for i = 2:nx-1;
      UNew(i) = (U(i+1) + U(i-1)) /2;
    end

    delta = max(abs(UNew - U), [], 'all');
    U = UNew;

    if(delta <= 1e-5)
      converged = true;
      break;
    end
  end

  residual = max(abs(U(1:end-2) - 2*U(2:end-1) + U(3:end)));

  if(~converged)
    warning('jacobi: достигнут лимит итераций (%d), max|dT| = %g, невязка = %g', ...
            maxIterations, delta, residual);
  end
end

