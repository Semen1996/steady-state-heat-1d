function T = gauss(T)
    % T = gauss(T)
    %
    % Gaussian elimination specialised for the tridiagonal matrix
    % of the interior-node system (diagonals 1, -2, 1).
    % The full matrix is not stored: zeros outside the three diagonals
    % stay zero during elimination and the off-diagonal entries equal 1,
    % so only the main diagonal d and the right-hand side b are kept.
    % No pivoting: the matrix is diagonally dominant, |d(i)| >= 1.
    % In effect this is the Thomas algorithm (thomas.m) with a = c = 1.
    %
    % Input:
    %   T - row vector of length nx; T(1) and T(end) are the boundary
    %       temperatures TL and TR, interior values are ignored
    %
    % Output:
    %   T - solution: the same vector with the interior nodes filled in


    nx = length(T);
    n = nx - 2;

    TL = T(1);
    TR = T(end);

    % Main diagonal and right-hand side
    d = -2 * ones(n, 1);
    b = zeros(n, 1);

    % Boundary conditions
    b(1) = -TL;
    b(end) = -TR;

    % Forward elimination
    for i = 2:n
        m = 1 / d(i-1);
        d(i) = d(i) - m;
        b(i) = b(i) - m * b(i-1);
    end

    % Back substitution
    TInner = zeros(n, 1);
    TInner(n) = b(n) / d(n);

    for i = n-1:-1:1
        TInner(i) = (b(i) - TInner(i+1)) / d(i);
    end

    T(2:end-1) = TInner.';
end

