clear all;
%% Добавление папки с методами
addpath(fullfile(pwd, 'methods', 'direct-methods'));
addpath(fullfile(pwd, 'methods', 'iterative-methods'));

%% Физические параметры
L = 1; % длина пластинки, [м]
TL = 100; % Температура слева, [град]
TR = 20; % Температура справа, [град]

%% Рачетные параметры
nx = 50;
dx = L / (nx - 1);
x = 0:dx:L;


% Коэффициенты
c2 = TL;
c1 = ( TR - TL ) / (L);

T_exact = c1*x + c2;

%% Начальное приближение
T = zeros(1,nx);
T(1) = TL;
T(end) = TR;


%% Решения

% 1. Прямые
T_gausse = gausse(T);
T_thomas = thomas(T);

% 2. Итерационные
T_jacoby = jacobi(T);
T_zeidel = gausse_zeidel(T);
T_sor = sor(T);




%% Демонстрация решения
figure(1), clf, hold on
plot(x, T_exact, 'r', 'LineWidth', 2);
plot(x, T_jacoby, 'b', 'LineWidth', 2);
plot(x, T_zeidel, 'g', 'LineWidth', 2);
plot(x, T_sor, 'm', 'LineWidth', 2);
plot(x, T_gausse, 'k', 'LineWidth', 2);
plot(x, T_thomas, 'c', 'LineWidth', 2);
grid on;
xlabel('x, m');
ylabel('T, grad');
legend('exact solution', 'Jacobi method', 'Gausse-Zeidel method', 'SOR method', 'Gausse method', 'Thomas method');
title('Temperature');
hold off
