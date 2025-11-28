# README

For the 1D Burgers equation problem, examples based on the following high-precision numerical method are provided:

- DG-RK3
- FD-RK3-WENO5
- FV-Euler-Godunov
- FV-RK3-WENO5

## Note

This project implements five C++ time integration frameworks, all supporting Forward Euler, SSP-RK3, as well as optional pre- and post-processing.
Each framework offers different trade-offs in flexibility, complexity, and performance, so all are retained:

- Runtime:
  - `solver_virtual.hpp`: virtual functions
  - `solver_stdfunc.hpp`: `std::function` (type erasure)

- Compile-time:
  - `solver_crtp.hpp`: CRTP
  - `solver_deducing.hpp`: deducing this (C++23)
  - `solver_template.hpp`: template parameters

Runtime approaches are simpler but slower; compile-time approaches are faster but more complex.

## Remark

C++23 is required:

- `std::format` (optional)
- `deducing this` (required for some examples, skip if unsupported)
