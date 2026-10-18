clear; close all;
%% Добавление папки с методами
addpath(fullfile(pwd, 'methods', 'direct-methods'));
addpath(fullfile(pwd, 'methods', 'iterative-methods'));

%% Физические параметры
L = 1; % длина стержня, [м]
TL = 100; % Температура слева, [град]
TR = 20; % Температура справа, [град]

%% Расчётные параметры
nx = 50;
dx = L / (nx - 1);
x = 0:dx:L;

% Параметры итерационных методов
tol = 1e-5; % порог max|T^(k+1) - T^(k)|
maxIterations = 10000;
omega = 1.5; % параметр релаксации SOR

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
T_jacoby = jacobi(T, tol, maxIterations);
T_zeidel = gausse_zeidel(T, tol, maxIterations);
T_sor = sor(T, tol, maxIterations, omega);




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
