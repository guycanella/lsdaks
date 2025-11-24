---
project: LSDA-Hubbard-Fortran
summary: Local Spin Density Approximation for the 1D Hubbard Model
author: Guilherme Canella
author_description: Computational Physics Researcher
email: guycanella@gmail.com
github: https://github.com/gcanella
project_github: https://github.com/gcanella/lsdaks
project_download: https://github.com/gcanella/lsdaks/releases
license: mit
docmark: <
predocmark: >
display: public
         protected
         private
source: true
graph: true
search: true
macro: TEST
       LOGIC=.true.
extra_mods: json_module: http://jacobwilliams.github.io/json-fortran/
            futility: http://cmacmackin.github.io
exclude_dir: ./test
             ./build
extra_filetypes: sh #

---

[TOC]

# LSDA-Hubbard-Fortran

A modern Fortran implementation of Local Spin Density Approximation (LSDA)
for solving the 1D Hubbard model using the Bethe Ansatz.

## Overview

This code implements density functional theory (DFT) in the local spin density approximation
for the one-dimensional Hubbard model. The exchange-correlation functional is computed
from exact Bethe Ansatz solutions using bicubic spline interpolation.

### Key Features

- **Exact Bethe Ansatz solver** using Newton-Raphson with analytical Jacobian
- **Exchange-correlation functional** via 2D bicubic spline interpolation
- **Self-consistent Kohn-Sham solver** with adaptive mixing
- **7 types of external potentials**: uniform, harmonic, barriers, disorder, quasiperiodic, impurities
- **High-performance linear algebra** using LAPACK (DSYEVD/ZHEEVD)
- **Comprehensive test suite** with 252+ tests (100% pass rate)
- **Production-ready**: Validated against C++ reference with energy agreement < 1e-8

## Physical Model

The 1D Hubbard Hamiltonian:

$$
H = -t \sum_{i,\sigma} (c^\dagger_{i\sigma} c_{i+1,\sigma} + h.c.)
    + U \sum_i n_{i\uparrow} n_{i\downarrow}
    + \sum_{i\sigma} V^{\text{ext}}_i n_{i\sigma}
$$

where:
- $t = 1$ (hopping parameter, energy unit)
- $U$ = Hubbard interaction (repulsive $U > 0$, attractive $U < 0$)
- $V^{\text{ext}}$ = External potential

## Code Architecture

### Module Hierarchy

```
types/              - Core data structures and constants
  ├── lsda_types      - System parameters, results types
  ├── lsda_constants  - Physical/numerical constants
  └── lsda_errors     - Error codes and handling

bethe_ansatz/       - Exact Bethe Ansatz solver
  ├── bethe_equations    - Lieb-Wu equations
  ├── nonlinear_solvers  - Newton-Raphson, Broyden
  ├── continuation       - Continuation in U parameter
  ├── table_io           - XC table I/O (ASCII/binary)
  └── bethe_tables       - Generate XC tables

xc_functional/      - Exchange-correlation functional
  ├── spline2d        - 2D bicubic spline interpolation
  └── xc_lsda         - LSDA functional interface

potentials/         - External potentials
  ├── potential_uniform      - Constant potential
  ├── potential_harmonic     - Harmonic trap
  ├── potential_impurity     - Impurity potentials
  ├── potential_random       - Disorder potentials
  ├── potential_barrier      - Rectangular barriers
  ├── potential_quasiperiodic - Aubry-André-Harper
  └── potential_factory      - Factory pattern

hamiltonian/        - Hamiltonian construction
  ├── hamiltonian_builder    - Tight-binding Hamiltonian
  └── boundary_conditions    - Open/periodic/twisted BC

diagonalization/    - Eigensolvers
  ├── lapack_wrapper       - LAPACK interface
  └── degeneracy_handler   - Handle degenerate states

density/            - Density calculation
  └── density_calculator   - Density from eigenstates

convergence/        - SCF convergence
  ├── convergence_monitor  - Track convergence
  ├── mixing_schemes       - Linear mixing
  └── adaptive_mixing      - Adaptive mixing parameter

kohn_sham/          - Main SCF loop
  └── kohn_sham_cycle     - Self-consistent solver

io/                 - Input/Output
  ├── input_parser    - Parse namelist input
  └── output_writer   - Write results
```

## Building and Testing

See the main [README.md](|page|/index.html) for build instructions.

## Physics Background

### Bethe Ansatz

The Bethe Ansatz provides exact solutions for the 1D Hubbard model via the Lieb-Wu equations:

$$
k_j L + \sum_{\alpha=1}^M \theta(k_j - \Lambda_\alpha) = 2\pi I_j
$$

$$
\sum_{j=1}^N \theta(k_j - \Lambda_\alpha) = 2\pi J_\alpha + \sum_{\beta=1}^M \Theta(\Lambda_\alpha - \Lambda_\beta)
$$

### DFT-LSDA Mapping

The Kohn-Sham equations are solved self-consistently:

$$
\left[-\nabla^2 + V_{\text{eff},\sigma}(r)\right] \phi_{i\sigma}(r) = \epsilon_{i\sigma} \phi_{i\sigma}(r)
$$

where $V_{\text{eff},\sigma} = V_{\text{ext}} + U n_{\bar{\sigma}} + V_{xc,\sigma}[n_\uparrow, n_\downarrow]$

## References

- E.H. Lieb and F.Y. Wu, *Phys. Rev. Lett.* **20**, 1445 (1968) - Original Bethe Ansatz
- F.H.L. Essler et al., *The One-Dimensional Hubbard Model* (Cambridge, 2005)
- K. Capelle and V.L. Campo, *Phys. Rep.* **528**, 91 (2013) - DFT for model Hamiltonians

## License

MIT License - see LICENSE file for details.

---

*Documentation generated with [FORD](https://github.com/Fortran-FOSS-Programmers/ford)*
