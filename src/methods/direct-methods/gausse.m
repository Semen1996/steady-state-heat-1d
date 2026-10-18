
function T = gausse(T)
    nx = length(T);
    n = nx - 2;

    TL = T(1);
    TR = T(end);

    % Главная диагональ и правая часть
    d = -2 * ones(n, 1);
    b = zeros(n, 1);

    % Учет граничных условий
    b(1) = -TL;
    b(end) = -TR;

    % Прямой ход метода Гаусса
    for i = 2:n
        m = 1 / d(i-1);
        d(i) = d(i) - m;
        b(i) = b(i) - m * b(i-1);
    end

    % Обратная подстановка
    Tint = zeros(n, 1);
    Tint(n) = b(n) / d(n);

    for i = n-1:-1:1
        Tint(i) = (b(i) - Tint(i+1)) / d(i);
    end

    T(2:end-1) = Tint.';
end

