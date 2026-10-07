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
graph_maxnodes: 40
graph_maxdepth: 3
search: true
macro: TEST
       LOGIC=.true.
extra_mods: json_module: http://jacobwilliams.github.io/json-fortran/
            futility: http://cmacmackin.github.io
src_dir: ./src
         ./app
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
- **Six potential-generator families** (ten selectable variants): uniform, harmonic, barriers, disorder, quasiperiodic, impurities
- **High-performance linear algebra** using LAPACK (DSTEVR/ZHEEVR)
- **Comprehensive test suite**: 20 Fortuno suites (current counts in the README, "Running Tests")
- **Validation**: analytic anchors and a historical C++ comparison, documented case by case in the README

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
  ├── nonlinear_solvers  - Newton-Raphson with analytical Jacobian
  ├── continuation       - Continuation in U parameter
  ├── lieb_wu_integral   - Thermodynamic-limit Lieb-Wu integral equations
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
  ├── potential_seed         - Seed handling for random potentials
  └── potential_factory      - Factory pattern

hamiltonian/        - Hamiltonian construction
  ├── hamiltonian_builder    - Tight-binding Hamiltonian
  └── boundary_conditions    - Open/periodic/twisted BC

diagonalization/    - Eigensolvers
  └── lapack_wrapper       - LAPACK interface (DSTEVR/DSYEVR/ZHEEVR)

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

Regenerate this documentation with `ford ford.md` from the repository root. The
graphs are bounded on purpose (`graph_maxnodes: 40`, `graph_maxdepth: 3` in the
header above): with `app/` in `src_dir`, unbounded graphs make FORD 7.0.12 hang
indefinitely at the graph stage, because the programs' call trees span the whole
project. Bounded, the run takes about 100 s. The header must contain only
`key: value` lines - a `#` line is not a comment there and silently ends it.

See the main [README.md](|page|/index.html) for build instructions.

## Physics Background

### Bethe Ansatz

The Bethe Ansatz provides exact solutions for the 1D Hubbard model via the Lieb-Wu equations:

$$
k_j L = 2\pi I_j - \sum_{\alpha=1}^M \theta(\sin k_j - \Lambda_\alpha)
$$

$$
\sum_{j=1}^N \theta(\Lambda_\alpha - \sin k_j) = 2\pi J_\alpha + \sum_{\beta \neq \alpha} \Theta(\Lambda_\alpha - \Lambda_\beta)
$$

with $\theta(x) = 2\arctan(x/u)$, $\Theta(x) = 2\arctan(x/2u)$ and $u = U/4$, the
convention of `bethe_equations.f90` (`theta(sin(k(j)) - Lambda(alpha), U)`).

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
