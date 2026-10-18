clear; close all;
%% Добавление папки с методами
addpath(fullfile(pwd, 'methods', 'direct-methods'));
addpath(fullfile(pwd, 'methods', 'iterative-methods'));

%% Физические параметры
L = 1; % длина стержня, [м]
TL = 100; % Температура слева, [°C]
TR = 20; % Температура справа, [°C]

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

TExact = c1*x + c2;

%% Начальное приближение
T = zeros(1,nx);
T(1) = TL;
T(end) = TR;


%% Решения

% 1. Прямые
TGauss = gauss(T);
TThomas = thomas(T);

% 2. Итерационные
TJacobi = jacobi(T, tol, maxIterations);
TSeidel = gauss_seidel(T, tol, maxIterations);
TSor = sor(T, tol, maxIterations, omega);




%% Демонстрация решения
figure(1), clf, hold on
plot(x, TExact, 'r', 'LineWidth', 2);
plot(x, TJacobi, 'b', 'LineWidth', 2);
plot(x, TSeidel, 'g', 'LineWidth', 2);
plot(x, TSor, 'm', 'LineWidth', 2);
plot(x, TGauss, 'k', 'LineWidth', 2);
plot(x, TThomas, 'c', 'LineWidth', 2);
grid on;
xlabel('x, m');
ylabel('T, ^\circC');
legend('exact solution', 'Jacobi method', 'Gauss-Seidel method', 'SOR method', 'Gauss method', 'Thomas method');
title('Temperature');
hold off
