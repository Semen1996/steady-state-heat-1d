# Steady-State Heat Conduction in a Rod (1D)

A small GNU Octave project that solves the one-dimensional steady-state heat equation (Laplace equation)

$$\frac{d^2 T}{d x^2} = 0, \qquad T(0) = T_L, \quad T(L) = T_R$$

with the finite difference method, and compares several ways of solving the resulting linear system against the exact linear profile.

The project has two parts that go together:

- **Article** — theory, derivation of the finite-difference scheme, and a walkthrough of every method with hand-computable examples, cost estimates and pitfalls:
  - [docs/physics-en.md](docs/physics-en.md) — English
  - [docs/physics.md](docs/physics.md) — Russian
- **Code** — [src/](src/) implements exactly the methods and formulas described in the article.

## Methods

| Type      | Method                                   | File |
|-----------|------------------------------------------|------|
| Direct    | Gaussian elimination (tridiagonal form)  | [gauss.m](src/methods/direct-methods/gauss.m) |
| Direct    | Thomas algorithm (TDMA)                  | [thomas.m](src/methods/direct-methods/thomas.m) |
| Iterative | Jacobi                                   | [jacobi.m](src/methods/iterative-methods/jacobi.m) |
| Iterative | Gauss–Seidel                             | [gauss_seidel.m](src/methods/iterative-methods/gauss_seidel.m) |
| Iterative | Successive over-relaxation (SOR)         | [sor.m](src/methods/iterative-methods/sor.m) |

## Running

Requires [GNU Octave](https://octave.org). Run the script from the `src/` directory (solver paths are added relative to the current folder):

```sh
cd src
octave --no-gui --eval "main"
```

Or open Octave, `cd` into `src/` and type `main`.

The script plots all numerical solutions on top of the exact solution `TExact`. To check a single method, compare it with the exact profile, e.g. `max(abs(TJacobi - TExact))`.

## Parameters

All parameters are set at the top of [src/main.m](src/main.m):

| Parameter       | Default | Meaning |
|-----------------|---------|---------|
| `L`             | `1`     | rod length, m |
| `TL`, `TR`      | `100`, `20` | boundary temperatures, °C |
| `nx`            | `50`    | number of grid nodes |
| `tol`           | `1e-5`  | stopping threshold for `max\|T^(k+1) - T^(k)\|` |
| `maxIterations` | `10000` | iteration limit |
| `omega`         | `1.5`   | SOR relaxation parameter, `0 < omega < 2` |

## Solver interface

Every solver takes a row vector `T` of length `nx`, reads the boundary temperatures from `T(1)` and `T(end)`, and returns the same vector with the interior nodes filled in:

```matlab
T = gauss(T);
T = thomas(T);

[T, iter, residual] = jacobi(T, tol, maxIterations);
[T, iter, residual] = gauss_seidel(T, tol, maxIterations);
[T, iter, residual] = sor(T, tol, maxIterations, omega);
```

- Missing or empty (`[]`) arguments fall back to the defaults: `tol = 1e-5`, `maxIterations = 10000`, `omega = 1.5`.
- `iter` is the number of iterations performed; `residual` is `max|T(i-1) - 2T(i) + T(i+1)|` over the interior nodes.
- Iterative methods stop when `max|T^(k+1) - T^(k)| <= tol` and issue a `warning` if this is not reached within `maxIterations`.

Note that the system matrix and right-hand side are not passed in: each solver has the Dirichlet-problem system (diagonals `1, -2, 1`) built in. Other boundary conditions (Neumann, Robin) or a heat source would require changing the solvers.

## Project structure

```
.
├── docs/
│   ├── physics.md        # article (Russian)
│   ├── physics-en.md     # article (English)
│   └── grid-1D.png
└── src/
    ├── main.m            # entry point: parameters, solving, plotting
    └── methods/
        ├── direct-methods/      # gauss.m, thomas.m
        └── iterative-methods/   # jacobi.m, gauss_seidel.m, sor.m
```
