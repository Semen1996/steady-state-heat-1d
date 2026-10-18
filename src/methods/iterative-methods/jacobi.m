function U = jacobi(U)
  nx = length(U);
  maxIterations = 10000;


  UNew = U;

  for n = 0:maxIterations
    for i = 2:nx-1;
      UNew(i) = (U(i+1) + U(i-1)) /2;
    end

    if(max(abs(UNew - U), [], 'all') <= 1e-5)
      U = UNew;
      break;
    end

    U = UNew;
  end
end

