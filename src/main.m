% Steady-state heat conduction in a rod: d²T/dx² = 0, T(0) = TL, T(L) = TR.
% The finite-difference linear system is solved by direct and iterative methods
% and the results are compared with the exact solution.

clear; close all;
%% Add the solver folders to the path
addpath(fullfile(pwd, 'methods', 'direct-methods'));
addpath(fullfile(pwd, 'methods', 'iterative-methods'));

%% Physical parameters
L = 1; % rod length, [m]
TL = 100; % left temperature, [°C]
TR = 20; % right temperature, [°C]

%% Grid parameters
nx = 50;
dx = L / (nx - 1);
x = 0:dx:L;

% Iterative method parameters
tol = 1e-5; % threshold for max|T^(k+1) - T^(k)|
maxIterations = 10000;
omega = 1.5; % SOR relaxation parameter


%% Initial guess
T = zeros(1,nx);
T(1) = TL;
T(end) = TR;


%% Solutions

% 0. Exact
TExact = ( TR - TL ) / (L)*x + TL;

% 1. Direct
TGauss = gauss(T);
TThomas = thomas(T);

% 2. Iterative
TJacobi = jacobi(T, tol, maxIterations);
TSeidel = gauss_seidel(T, tol, maxIterations);
TSor = sor(T, tol, maxIterations, omega);




%% Plot the solutions
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
