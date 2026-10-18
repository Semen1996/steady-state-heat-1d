function [T, iter, residual] = jacobi(T, tol, maxIterations)
  % [T, iter, residual] = jacobi(T, tol, maxIterations)
  %
  % Метод Якоби для системы T(i-1) - 2T(i) + T(i+1) = 0, i = 2..nx-1
  % Новое значение в узле считается по значениям
  % соседей с предыдущей итерации.
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


  nx = length(T);
  if(nargin < 2 || isempty(tol))
    tol = 1e-5;
  end
  if(nargin < 3 || isempty(maxIterations))
    maxIterations = 10000;
  end


  TNew = T;
  converged = false;

  for iter = 1:maxIterations
    for i = 2:nx-1;
      TNew(i) = (T(i+1) + T(i-1)) /2;
    end

    delta = max(abs(TNew - T), [], 'all');
    T = TNew;

    if(delta <= tol)
      converged = true;
      break;
    end
  end

  residual = max(abs(T(1:end-2) - 2*T(2:end-1) + T(3:end)));

  if(~converged)
    warning('jacobi: достигнут лимит итераций (%d), max|dT| = %g, невязка = %g', ...
            maxIterations, delta, residual);
  end
end

