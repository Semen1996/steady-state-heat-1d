function U = gausse_zeidel(U)
  maxIterations = 10000;
  nx = length(U);

  UOld = U;

  for n = 0:maxIterations
    for i = 2:nx-1;
      U(i) = (U(i+1) + U(i-1)) /2;
    end

    if(max(abs(UOld - U), [], 'all') <= 1e-5)
      break;
    end

    UOld = U;
  end


end

