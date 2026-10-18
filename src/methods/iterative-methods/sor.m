function U = sor(U)
  nx = length(U);
  maxIterations = 10000;
  omega = 1.5;

  UOld = U;

  for n = 0:maxIterations
    for i = 2:nx-1;
      UStar = (U(i+1) + U(i-1)) / 2;
      U(i) = U(i) + omega*(UStar - U(i));
    end

    if(max(abs(UOld - U), [], 'all') <= 1e-5)
      break;
    end

    UOld = U;
  end
end

