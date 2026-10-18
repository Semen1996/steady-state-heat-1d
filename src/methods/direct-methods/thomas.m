function T = thomas(T)
    % T = thomas(T)
    %
    % Thomas algorithm (tridiagonal matrix algorithm) for the system
    % T(i-1) - 2T(i) + T(i+1) = 0, i = 2..nx-1.
    % Only the three diagonals a, b, c and the right-hand side d are stored;
    % the forward sweep eliminates the lower diagonal, the backward sweep
    % recovers the solution.
    %
    % Input:
    %   T - row vector of length nx; T(1) and T(end) are the boundary
    %       temperatures TL and TR, interior values are ignored
    %
    % Output:
    %   T - solution: the same vector with the interior nodes filled in


    nx = length(T);
    n = nx - 2;

    % Boundary temperatures
    TL = T(1);
    TR = T(end);

    % Coefficients of the tridiagonal system
    a = ones(n, 1);       % Lower diagonal
    b = -2 * ones(n, 1);  % Main diagonal
    c = ones(n, 1);       % Upper diagonal
    d = zeros(n, 1);      % Right-hand side

    % Boundary conditions
    a(1) = 0;
    c(end) = 0;
    d(1) = -TL;
    d(end) = -TR;

    % Forward sweep
    for i = 2:n
        m = a(i) / b(i-1);
        b(i) = b(i) - m * c(i-1);
        d(i) = d(i) - m * d(i-1);
    end

    % Backward sweep
    TInner = zeros(n, 1);
    TInner(n) = d(n) / b(n);

    for i = n-1:-1:1
        TInner(i) = (d(i) - c(i) * TInner(i+1)) / b(i);
    end

    % Store the interior temperatures
    T(2:end-1) = TInner.';
end

