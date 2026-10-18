function [T, iter, residual] = gauss_seidel(T, tol, maxIterations)
  % [T, iter, residual] = gauss_seidel(T, tol, maxIterations)
  %
  % Метод Гаусса-Зейделя для системы T(i-1) - 2T(i) + T(i+1) = 0,
  % i = 2..nx-1. Узлы обходятся слева направо,
  % уже пересчитанное значение левого соседа используется сразу.
  %
  % Вход:
  %   T             - вектор-строка длины nx: начальное приближение;
  %                   T(1) и T(end) - граничные температуры TL и TR
  %   tol           - порог max|T^(k+1) - T^(k)| (по умолчанию 1e-5)
  %   maxIterations - наибольшее число итераций (по умолчанию 10000)
  %
  % Выход:
  %   T        - решение (граничные значения не меняются)
  %   iter     - число выполненных итераций
  %   residual - невязка max|T(i-1) - 2T(i) + T(i+1)| по внутренним узлам
  %
  % Если за maxIterations итераций критерий max|T^(k+1) - T^(k)| <= tol
  % не выполнен, выдаётся предупреждение.


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
    warning('gauss_seidel: достигнут лимит итераций (%d), max|dT| = %g, невязка = %g', ...
            maxIterations, delta, residual);
  end
end

