# LSDA-Hubbard-Fortran

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Fortran](https://img.shields.io/badge/Fortran-2008%2F2018-734f96?logo=fortran)](https://fortran-lang.org)
[![fpm](https://img.shields.io/badge/fpm-compatible-brightgreen)](https://fpm.fortran-lang.org)

A modern Fortran implementation of **Local Spin Density Approximation (LSDA)** for solving the one-dimensional Hubbard model using the Bethe Ansatz. The code is modular, covered by an extensive test suite, and validated against the C++ reference on the cases listed under [Validation](#validation) — which is a narrower claim than "validated in general": several features carry documented caveats, flagged with ⚠️ where they appear.

## Table of Contents

- [Features](#features)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Directory Structure](#directory-structure)
- [External Potentials](#external-potentials)
- [Input Format](#input-format)
- [Output Files](#output-files)
- [Running Tests](#running-tests)
- [Documentation](#documentation)
- [Physical Background](#physical-background)
- [Performance](#performance)
- [Validation](#validation)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)

## Features

- ✅ **Exact Bethe Ansatz solver** using Newton-Raphson with analytical Jacobian
- ✅ **Exchange-correlation functional** interpolated with bicubic splines in density and magnetization on pre-computed tables
- ✅ **Self-consistent Kohn-Sham solver** with adaptive mixing for stability
- ✅ **Six potential-generator families** (ten selectable variants): uniform, harmonic, impurities, disorder, barriers, and quasiperiodic modulation
- ✅ **High-performance linear algebra** using LAPACK (DSTEVR/ZHEEVR)
- ✅ **Three boundary conditions**: open, periodic, twisted
- ✅ **20 test suites** with **1,455 `call check` assertions** (`grep -c "call check(" test/*.f90`)
- ✅ **Measured C++ comparisons** and deliberate compatibility differences documented below

## System Requirements

### Required

- **Fortran compiler** with Fortran 2008/2018 support:
  - GCC `gfortran` ≥ 9.0
  - Intel `ifort` ≥ 19.0
  - LLVM `flang` (recent versions)

- **LAPACK/BLAS libraries**:
  - Linux: `liblapack-dev`, `libblas-dev`
  - macOS: Included in Accelerate framework (automatic)
  - Windows: Intel MKL or OpenBLAS

- **Fortran Package Manager (fpm)**:
  ```bash
  # Install via conda (recommended)
  conda install -c conda-forge fpm

  # Or via pip
  pip install fpm

  # Or download binary from https://github.com/fortran-lang/fpm/releases
  ```

### Optional

- **Python 3.7+** (for analysis scripts and benchmark tools)
  - `numpy`, `matplotlib` (for plotting)

- **FORD** (for generating HTML documentation):
  ```bash
  pip install ford
  ```

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/yourusername/lsdaks.git
cd lsdaks
```

### 2. Install System Dependencies

#### Ubuntu/Debian
```bash
sudo apt-get update
sudo apt-get install gfortran liblapack-dev libblas-dev
```

#### macOS (with Homebrew)
```bash
brew install gcc  # Includes gfortran
# LAPACK/BLAS included in Accelerate framework (no action needed)
```

#### macOS (with MacPorts)
```bash
sudo port install gcc12  # Or gcc13, gcc14
sudo port select --set gcc mp-gcc12
```

### 3. Install Fortran Package Manager

```bash
# Via conda (recommended)
conda install -c conda-forge fpm

# Or manually download from
# https://github.com/fortran-lang/fpm/releases
```

### 4. Build the Project

```bash
# Development build
fpm build

# Optimized build for production
fpm build --profile release --flag "-O3 -march=native"
```

**Note**: The first build will download and compile dependencies (Fortuno testing framework). This may take a few minutes.

## Quick Start

### Running a Simple Calculation

1. Create an input file `input.txt`:

```fortran
&system
  L = 100              ! Number of lattice sites
  Nup = 50             ! Number of spin-up electrons
  Ndown = 50           ! Number of spin-down electrons
  U = 4.0              ! Hubbard interaction
  bc = 'open'          ! Boundary condition: 'open', 'periodic', or 'twisted'
/

&potential
  potential_type = 'uniform'  ! Uniform potential
  V0 = 0.0                    ! Potential value
/

&scf
  max_iter = 1000             ! Maximum SCF iterations
  density_tol = 1.0e-6        ! Density convergence tolerance
  mixing_alpha = 0.05         ! Mixing parameter (0 < α ≤ 1)
  use_adaptive_mixing = .true. ! Use adaptive mixing
  verbose = .false.           ! Print iteration details
/

&output
  output_prefix = 'results'
  save_density = .true.
  save_eigenvalues = .true.
  save_wavefunction = .false.
/
```

2. Run the calculation:

```bash
fpm run --profile release --flag "-O3 -march=native" lsdaks -- --input input.txt
```

3. Results are saved to:
   - `results_summary.txt` - Final energy and convergence info
   - `results_density.dat` - Site-resolved densities
   - `results_eigenvalues.dat` - Kohn-Sham eigenvalues
   - `results_convergence.dat` - SCF convergence history

### Using Example Inputs

The `examples/` directory contains pre-configured input files:

```bash
# Minimal example
fpm run lsdaks -- --input examples/input_minimal.txt

# Half-filling (n=1.0) validation case
fpm run lsdaks -- --input examples/input_halffilling.txt

# Harmonic trap
fpm run lsdaks -- --input examples/input_harmonic_trap.txt

# Strong coupling regime
fpm run lsdaks -- --input examples/input_strong_coupling.txt

# Twisted boundary conditions
fpm run lsdaks -- --input examples/input_twisted_bc.txt
```

## Directory Structure

```
lsdaks/
├── app/                          # Executable programs
│   ├── main.f90                  # Main LSDA solver
│   ├── convert_tables.f90        # XC table format converter
│   └── generate_table.f90        # Generate thermodynamic-limit XC tables via Bethe Ansatz
│
├── src/                          # Source code
│   ├── types/                    # Core data structures
│   │   ├── lsda_types.f90        # System parameters, results types
│   │   ├── lsda_constants.f90    # Physical/numerical constants
│   │   └── lsda_errors.f90       # Error handling
│   │
│   ├── bethe_ansatz/             # Bethe Ansatz solver
│   │   ├── bethe_equations.f90   # Lieb-Wu equations
│   │   ├── nonlinear_solvers.f90 # Newton-Raphson
│   │   ├── continuation.f90      # Continuation in U
│   │   ├── table_io.f90          # Table I/O (ASCII/binary)
│   │   └── bethe_tables.f90      # Generate XC tables
│   │
│   ├── xc_functional/            # Exchange-correlation
│   │   ├── spline2d.f90          # bicubic XC interpolation
│   │   └── xc_lsda.f90           # LSDA functional interface
│   │
│   ├── potentials/               # External potentials
│   │   ├── potential_uniform.f90
│   │   ├── potential_harmonic.f90
│   │   ├── potential_impurity.f90
│   │   ├── potential_random.f90
│   │   ├── potential_barrier.f90
│   │   ├── potential_quasiperiodic.f90
│   │   └── potential_factory.f90
│   │
│   ├── hamiltonian/              # Hamiltonian construction
│   │   ├── hamiltonian_builder.f90
│   │   └── boundary_conditions.f90
│   │
│   ├── diagonalization/          # Eigensolvers
│   │   └── lapack_wrapper.f90
│   │
│   ├── density/                  # Density calculation
│   │   └── density_calculator.f90
│   │
│   ├── convergence/              # SCF convergence
│   │   ├── convergence_monitor.f90
│   │   ├── mixing_schemes.f90
│   │   └── adaptive_mixing.f90
│   │
│   ├── kohn_sham/                # Main SCF loop
│   │   └── kohn_sham_cycle.f90
│   │
│   └── io/                       # Input/Output
│       ├── input_parser.f90
│       └── output_writer.f90
│
├── test/                         # Test suite (20 suites, 1,455 assertions)
│   ├── test_bethe_equations.f90
│   ├── test_nonlinear_solvers.f90
│   ├── test_continuation.f90
│   ├── test_spline2d.f90
│   ├── test_xc_lsda.f90
│   ├── test_potentials.f90
│   ├── test_hamiltonian_builder.f90
│   ├── test_lapack_wrapper.f90
│   ├── test_density_calculator.f90
│   ├── test_convergence_monitor.f90
│   ├── test_kohn_sham_cycle.f90
│   └── ...
│
├── data/                         # Data files
│   ├── tables/                   # XC functional tables
│   │   └── fortran_native/       # Binary format (fast loading)
│
├── examples/                     # Example input files
│   ├── input_minimal.txt
│   ├── input_halffilling.txt
│   ├── input_harmonic_trap.txt
│   ├── input_strong_coupling.txt
│   └── input_twisted_bc.txt
│
├── benchmark_results/            # Validation data
│   ├── benchmark_table.md
│   ├── detailed_report.txt
│   └── *.png                     # Comparison plots
│
├── scripts/                      # Utility scripts
├── fpm.toml                      # FPM build configuration
├── ford.md                       # FORD documentation config
├── CLAUDE.md                     # AI assistant guidance
├── PROJECT_CONTEXT.md            # Detailed technical docs (Portuguese)
└── README.md                     # This file
```

### XC Tables and the Table Generator

The SCF reads pre-computed exchange-correlation tables from `data/tables/fortran_native/`
(25 values of U, converted from the C++ reference; Hartree already subtracted). These are the
only tables validated for production runs.

`generate_xc_table` solves the thermodynamic-limit Lieb-Wu integral equations on its
configured, graded `(n, m)` grid and writes the native table format consumed by the SCF. The
default grid has 75 density rows and 202 magnetization nodes per row, graded on **both** axes
(`density_grid` and `magnetization_grid` in `src/bethe_ansatz/bethe_tables.f90`); it begins
at `n_min = 0.02` and therefore does not cover the entire physical density--magnetization
triangle. Its parameters and quadrature orders can be adjusted on the command line. A U=4
table on this configured grid can be generated with the default settings. The executable
refuses to write a table containing NaN or Inf, reports the offending grid points, and exits
with status 1.

**How far the generated tables are validated.** The generator was compared against the
converted C++ reference tables on the **15 integer values of U**: worst `|Δexc| = 8.5e-7`,
zero nodes outside 1e-6. For U=4, along the path the SCF actually consumes
(`xc_lsda_init` + `get_exc`/`get_vxc`, 10099 nodes, corner n=1,m=1 excluded), worst
`|Δexc| = 2.75e-7` and worst `|ΔVxc| = 7.67e-5`, zero nodes outside tolerance. No
equivalence is claimed for non-integer U: the 10 reference tables in that range are
internally inconsistent and were rejected as an oracle (see the "Resultado" of phase 4.5 in
`NEXT_STEPS_REPORT.md`).

The default output directory is the SCF table directory. To protect its existing reference
tables, the generator refuses to overwrite an existing U table unless `--force` is supplied;
use `--output <directory>` when creating a separate table set.

**An XC table is mandatory for every SCF run.** The supported range is
`1 <= |U| <= 20`: the lower end is `bethe_tables::U_TABLE_MIN`
(`src/bethe_ansatz/bethe_tables.f90:153`), the upper end is simply the largest shipped
table (`data/tables/fortran_native/xc_table_u20.00.dat`) — nothing in the code rejects
`|U| > 20`, but there is no table there, so `U = 25` fails with the same "XC table not
found!" error as `U = 0`.

The XC *physics* of the non-interacting limit is already in the code: both `get_exc` and
`get_vxc` short-circuit to zero for `|U| < U_SMALL = 1e-9`
(`src/xc_functional/xc_lsda.f90:283` and `:473`). What blocks `U = 0` is only the table
**loading** in `xc_lsda_init`, which the SCF performs unconditionally: `app/main.f90:94`
resolves the table for `|U|` and, when the file is absent, `app/main.f90:97-110` prints
"ERROR: XC table not found!" and stops with status 1. Table generation is refused for
`U = 0` as well (`src/bethe_ansatz/bethe_tables.f90:537`), so `U = 0` cannot be run today.
The failure mode is also cosmetically confusing: the file name is
built with the `F0.2` descriptor (`app/main.f90:629`), which writes `0.0` as `.00` without
the leading zero, so the executable looks for `xc_table_u.00.dat`, and the command the error
message suggests (`fpm run generate_xc_table -- --U .00`) is itself refused by the
generator. `F0.2` produces the correct name for every `U >= 1`, so the formatting defect only
shows up in the range that is already unsupported. A future task only needs to let
`xc_lsda_init` skip the file read when `|U| < U_SMALL` — the zero-XC branches downstream
already exist; that change is outside T17.

## External Potentials

The code exposes 10 selectable variants from six potential-generator families:

### 1. Uniform Potential
```fortran
&potential
  potential_type = 'uniform'
  V0 = 0.0  ! Constant value
/
```
- **Formula**: `V(i) = V₀`
- **Use case**: Homogeneous systems, baseline for testing

### 2. Harmonic Trap
```fortran
&potential
  potential_type = 'harmonic'
  spring_constant = 0.02  ! Trap strength k
/
```
- **Formula**: `V(i) = k × (i - i_center)²`
- **Center**: `i_center = (L+1)/2` (middle of chain)
- **Use case**: Optical traps in cold atoms, shell structure

### 3. Single Impurity
```fortran
&potential
  potential_type = 'impurity_single'
  V0 = 2.0          ! Impurity strength
  pot_center = 50.0 ! Impurity site (1-indexed)
/
```
- **Formula**: `V(i) = V_imp` if `i = i_pos`, else `V(i) = 0`
- **Use case**: Point defects, Kondo physics

### 4. Multiple Impurities
```fortran
&potential
  potential_type = 'impurity_multiple'
  V0 = 2.0
  imp_positions_str = '25, 50, 75'
/
```
- **Use case**: Multiple defects, disorder modeling

### 5. Random Impurities (by Concentration)
```fortran
&potential
  potential_type = 'impurity'
  V0 = 2.0              ! Impurity strength
  concentration = 10.0  ! 10% of sites have impurities
  pot_seed = 12345      ! For reproducibility (-1 = random)
/
```
- **Formula**: Randomly places `N_imp = round(concentration × L / 100)` impurities
- **Example**: `L=100, concentration=10.0` → 10 impurities at random positions
- **Use case**: Dilute disorder, Anderson localization

### 6. Uniform Disorder
```fortran
&potential
  potential_type = 'random_uniform'
  disorder_strength = 2.0  ! Width W
  pot_seed = 12345
/
```
- **Formula**: `V(i) ~ Uniform[-W, +W]`
- **Mean**: `⟨V⟩ = 0`, **Variance**: `σ² = W²/3`
- **Use case**: Box disorder, Anderson localization

### 7. Gaussian Disorder
```fortran
&potential
  potential_type = 'random_gaussian'
  disorder_strength = 1.0  ! Std deviation σ
  pot_seed = 12345
/
```
- **Formula**: `V(i) ~ Normal(0, σ²)`
- **Use case**: Thermal/quantum fluctuations

### 8. Single Barrier
```fortran
&potential
  potential_type = 'barrier_single'
  V0 = 5.0       ! Barrier height
  position = 50  ! Reference site (1-indexed)
  width = 4      ! Number of covered sites
/
```
- **Formula**: `V(i) = V_b` if `i_start ≤ i ≤ i_end`, else `V(i) = 0`
- **Placement**: an odd width is centred on `position`. For an even width,
  the covered interval is `[position - width/2, position + width/2 - 1]`;
  e.g. `position = 50`, `width = 4` covers sites 48–51, with `position` as
  the right-hand one of the two central sites.
- **Use case**: Quantum tunneling, scattering

### 9. Double Barrier (Quantum Well)
```fortran
&potential
  potential_type = 'barrier_double'
  V0 = 5.0            ! Barrier height
  barrier_width = 5.0
  well_depth = -3.0
  well_width = 20.0
/
```
- **Geometry**: `[Barrier] [Well] [Barrier]`
- **Use case**: Resonant tunneling, quasi-bound states

### 10. Quasiperiodic (Aubry-André-Harper)
```fortran
&potential
  potential_type = 'quasiperiodic'
  aah_lambda = 2.0  ! Potential strength λ
  aah_beta = 0.618  ! Frequency β (typically (√5-1)/2)
  aah_phi = 0.0     ! Phase φ, in radians
/
```
- **Formula**: `V(i) = λ cos(2π β i + φ)`
- **Use case**: Anderson localization transition, topological physics

## Input Format

Input files use Fortran namelists (case-insensitive and order-independent). The table below lists the keys used in practice; it is not exhaustive (for example `pot_width` is also accepted in `/potential`, see `src/io/input_parser.f90`). Names that belong to no group are rejected. `phase` is supplied in units of π and converted to radians internally.

| Group | Key | Type | Default | Applies to |
|---|---|---|---|---|
| system | `L` | integer | 10 | all |
| system | `Nup`, `Ndown` | integer | 5, 5 | all |
| system | `U` | real | 4.0 | all |
| system | `bc` | character | `periodic` | all (`open`, `periodic`, `twisted`) |
| system | `phase` | real | 0.0 | `twisted`; input unit π |
| system | `table_dir` | character | `data/tables/fortran_native` | all |
| potential | `potential_type` | character | `uniform` | all |
| potential | `V0` | real | 0.0 | uniform, impurities, barriers |
| potential | `spring_constant` | real | 0.001 | harmonic |
| potential | `pot_center` | real | 0.0 | impurity_single |
| potential | `imp_positions_str` | character | empty | impurity_multiple |
| potential | `concentration` | real | 50.0 | impurity (random placement) |
| potential | `pot_seed` | integer | -1 | impurity, random_uniform, random_gaussian |
| potential | `disorder_strength` | real | 2.0 | random_uniform, random_gaussian |
| potential | `position`, `width` | integer | 50, 5 | barrier_single |
| potential | `barrier_width`, `well_depth`, `well_width` | real | 3.0, -3.0, 20.0 | barrier_double |
| potential | `aah_lambda`, `aah_beta`, `aah_phi` | real | 1.0, 0.6180339887498948, 0.0 | quasiperiodic; φ is radians |
| potential | `position1`, `width1`, `position2`, `width2` | integer | 35, 3, 65, 3 | deprecated; do not use |
| scf | `max_iter` | integer | `ITER_MAX` | all |
| scf | `density_tol` | real | `SCF_DENSITY_TOL` | diagnostic only |
| scf | `energy_tol`, `potential_tol` | real | `SCF_ENERGY_TOL`, `SCF_POTENTIAL_TOL` | convergence |
| scf | `mixing_alpha` | real | `MIX_ALPHA` | all; new-potential weight |
| scf | `verbose`, `store_history`, `use_adaptive_mixing` | logical | true, true, true | all |
| scf | `xc_smoothing_width` | real | 0.0 | opt-in smoothing near n=1 |
| output | `output_prefix` | character | `lsda_output` | all |
| output | `save_density`, `save_eigenvalues`, `save_wavefunction` | logical | true, true, false | all |

The selectable values of `potential_type` are `uniform`, `harmonic`, `impurity`, `impurity_single`, `impurity_multiple`, `random_uniform`, `random_gaussian`, `barrier_single`, `barrier_double`, and `quasiperiodic`. These are the strings **accepted in the namelist**
(`app/main.f90`), which is not the same list as the ten identifiers **recognised by
`get_potential_info`** (`src/potentials/potential_factory.f90:160-178`, documented in
`PROJECT_CONTEXT.md`): the namelist accepts the legacy alias `impurity` and routes random
placement through it, while that inventory knows `impurity_random` instead. Of those ten,
`create_potential` (`:78-146`) builds eight directly and rejects `impurity_multiple` and
`impurity_random` with `ERROR_INVALID_INPUT`; those two are served by specialised routines
in the `app/main.f90` dispatch. `distribution`, `twisted_phase`, `output_file`, `write_density`, `write_eigenvalues`, and `write_convergence_history` are not namelist keys.

### Complete Example

```fortran
&system
  L = 100              ! Lattice sites
  Nup = 50             ! Spin-up electrons
  Ndown = 50           ! Spin-down electrons
  U = 4.0              ! Hubbard U
  bc = 'open'          ! Boundary: 'open', 'periodic', 'twisted'
  phase = 0.0          ! Input phase in units of π; solver converts to radians
/

&potential
  potential_type = 'harmonic'
  spring_constant = 0.02
/

&scf
  max_iter = 10000
  density_tol = 1.0e-6
  energy_tol = 1.0e-8
  mixing_alpha = 0.05          ! Linear mixing (0 < α ≤ 1)
  use_adaptive_mixing = .true. ! Adjust α dynamically
  verbose = .false.            ! Print each iteration
/

&output
  output_prefix = 'results'
  save_density = .true.
  save_eigenvalues = .true.
  save_wavefunction = .false.
/
```

### Notes
- **Mixing convention**: `α = 0.05` means 5% new, 95% old (conservative)
- **Twisted BC**: `phase` in units of π (e.g., `phase = 0.5` → π/2)
- **Adaptive mixing**: Automatically adjusts `α` when convergence stalls
- **XC smoothing** (`xc_smoothing_width = w` in `&scf`, default `0`): replaces
  the BALDA discontinuity of `V_xc` at `n = 1` by a linear ramp of half-width
  `w`. The default stays `0` (exact C++ functional) because `w > 0` changes the
  Kohn--Sham potential without smoothing `E_xc` itself: for `w > 0`,
  `V_xc != δE_xc/δn` by a finite amount and the reported energy is not
  stationary at the fixed point. Note that `w = 0` does **not** restore an exact
  variational principle either: `e_xc`, `V_xc^up` and `V_xc^dn` are three
  independent bicubic splines over three independently tabulated columns
  (`src/xc_functional/xc_lsda.f90:36-38`, initialised separately at `:228`,
  `:232` and `:235`), so `V_xc` is never the analytic derivative of the `e_xc`
  spline. With `w = 0` the fixed point is stationary only up to the internal
  inconsistency between the splines, `‖V_xc^σ − ∂e_xc/∂n_σ‖`, measured directly
  on the shipped `xc_table_u4.00.dat` splines (finite differences, `h = 1e-5`,
  swept over `(n, m)`): of order `5e-4` away from half filling
  (`max |∂e_xc/∂n↑ − V_xc^↑| = 4.66e-4` for `|n - 1| > 0.05`), and with **no
  small bound at `n = 1`**, where BALDA has a discontinuity that the spline
  smooths out (`8.28e-1` at `n = 1.00, m = 0.975`). Do not confuse this with the
  generator-vs-reference table agreement quoted above (`|Δexc| = 2.75e-7`,
  `|ΔVxc| = 7.67e-5`): that is the distance between two tables and says nothing
  about how far `V_xc` is from the derivative of the `e_xc` spline. The `n = 1`
  regime is precisely where the smoothing below is recommended, so the reported
  energy there should not be read as variationally stationary.
  Reference results must still be produced with
  `w = 0`, because that is the only setting that reproduces the C++
  functional. Use smoothing as an explicit,
  per-case opt-in for systems whose density sits on `n = 1` (Mott plateau in a
  trap, double-barrier well); it is echoed in the output header when active.
- **Near-degenerate Fermi levels** (deliberate divergence from the C++):
  two consecutive Kohn-Sham levels closer than `DEG_TOL = 1e-10` share the
  open-shell occupation equally, matching the C++ `update_degen` rule for an
  exactly degenerate block; but
  between `1e-10` and `DEG_TOL_UPPER = 1e-6` the sharing fades out with a C¹
  weight instead of switching off abruptly (`compute_occupations`). This is an
  effective smearing over that `~1e-6 t` interval, not a strictly sharp
  zero-temperature occupation at every finite gap. Particle number is
  conserved algebraically (up to floating-point rounding), `0 <= occ <= 1`,
  and this occupation-rule equivalence does not imply bit-for-bit equivalence
  of the SCF result.

  Two consequences worth stating explicitly:

  1. **Energy cost.** Inside the transition window the occupations are
     fractional, so `E_band = Σ occ·ε` sits above the strict Aufbau sum by up to
     roughly `O(g(g-1)·1e-6)` for a `g`-fold near-degenerate shell: the shell
     spans at most `(g-1)·DEG_TOL_UPPER`, but the redistributed charge is `O(g)`.
     For the common case `g = 2` this is about `1e-6` (in practice `<= 5e-7`).
     That is roughly three orders of magnitude above the `5e-10` resolution
     floor of the C++ comparison below, measured on cases with no level
     inside the window. **The two numbers must not be treated as competing
     accuracy estimates**: the `5e-10` floor is the limit of what the C++
     comparison can resolve where the occupation rule is inactive, the
     `1e-6` figure is the deliberate divergence where it is active.
  2. **It is not a standard `f(ε - μ)` smearing.** The weight of level `j` is a
     *product* of link weights along the chain of consecutive neighbours
     (`src/density/density_calculator.f90:167-176`), and the chain only stops at
     the first fully open link. A level three small gaps away from the Fermi
     level can therefore still be pulled into the shared pool, which no
     function of `ε - μ` alone would do.

  The upper edge is absolute: in very large PBC systems
  (roughly `L >= 1e4`), distinct levels
  can enter the transition interval. A tunnel doublet returned by LAPACK in a
  localised basis can therefore no longer flip a whole electron into one arm of
  a symmetric trap as the gap fluctuates around `1e-10`, so `n(i) = n(L+1-i)`
  is preserved.

## Output Files

After a successful run, the following files are created:

### 1. Summary File (`results_summary.txt`)
```
System Parameters:
  L (sites):        100
  N_up:             50
  N_down:           50
  N_total:          100
  U:                4.0000
  BC:               open

SCF Convergence:
  Status:           ✓ CONVERGED
  Iterations:       127
  Final |Δn|:       8.3421E-07
  Final Total Energy:    -45.234567890123
  Final Energy per site:   -0.452345678901

Density Check:
  ∫n_up dx:         50.000000
  ∫n_down dx:       50.000000
  ∫n_total dx:      100.000000
  Expected N:       100.000000
  Error:            2.8422E-14
```

### 2. Density Profile (`results_density.dat`)
```
# Columns: site  n_up  n_down  n_total
     1    4.8566E-01    4.8566E-01    9.7134E-01
     2    5.0894E-01    5.0894E-01    1.0179E+00
     3    5.1799E-01    5.1799E-01    1.0360E+00
   ...
```

### 3. Eigenvalues (`results_eigenvalues.dat`)
```
# Columns: index  spin  eigenvalue  occupied
     1      up   -3.1234567890E+00     yes
     2      up   -2.9876543210E+00     yes
   ...
    50      up   -0.5432109876E+00     yes
    51      up    0.1234567890E+00     no
   ...
```

**The file holds up to `L` records per spin, and the count varies.** The SCF
diagonalizes only the occupied levels plus a small buffer for the Fermi shell,
so the levels above that window are never computed and are not written. A file
with, say, 55 spin-up records for `L = 1000` is complete, not truncated: read
whatever records are present and treat the missing levels as *not computed*
rather than as missing data or a short write. The record count also varies with
the filling and, when a near-degenerate shell forces the window to grow, between
runs of the same system.

### 4. Convergence History (`results_convergence.dat`)
```
# Columns: iteration  energy  density_error  mixing_alpha
     1   -40.123456    5.6789E-02    0.0500
     2   -42.345678    3.4567E-02    0.0500
   ...
   127   -45.234567    8.3421E-07    0.0500
```

## Running Tests

The project has 20 explicitly registered suites and 1,455 assertions as of 2026-09-20 (`grep -c "call check(" test/*.f90`). Running `fpm test --profile release` on that date executed 351 Fortuno test cases across the 20 suites with 0 failures. Assertion counts are source inventory, not a claim about coverage; the per-suite `Total:` line that Fortuno prints is that suite's case count, not the project total.

> ⚠️ **Always run the tests with `--profile release`.**
>
> With gfortran 16, the default (debug) profile aborts 18 of the 20 suites inside
> the Fortuno test driver before a single test runs:
>
> ```
> At line 233 of file build/dependencies/fortuno/src/fortuno/testdriver.f90
> Fortran runtime error: Index '1' of dimension 1 of array 'this...%suiteresults'
> outside of expected range (0:0)
> ```
>
> This is debug-only `-fcheck=bounds` tripping over a zero-sized array inside
> the Fortuno dependency, not a defect in this project. Pinning Fortuno to its
> only published tag (`v0.1.0`) was tested and does **not** avoid the abort, so
> the dependency is left unpinned and the release profile is the supported way
> to run the suite until Fortuno fixes the zero-sized-array access.

### Run All Tests
```bash
fpm test --profile release
```

### Run Specific Test Suites
```bash
fpm test --profile release test_bethe_equations
fpm test --profile release test_kohn_sham_cycle
fpm test --profile release test_potentials
fpm test --profile release test_nonlinear_solvers
```

All test programs are listed explicitly as `[[test]]` blocks in `fpm.toml`
(`auto-tests` is disabled) so the test inventory is auditable. When adding a
test file under `test/`, add a matching `[[test]]` block.

### Test Coverage

The suite inventory is the 20 `[[test]]` blocks in `fpm.toml`. Run `fpm test --profile release` for the current result; do not infer a stable total from historical documentation.

### Validation Tests

Use `scripts/build_cpp_reference.sh` to build the reference into `build/cpp/`; no `run_all_tests.sh`, `run_cpp_tests.sh`, or `compare_energies.py` exists in this repository.

## Documentation

### Generate HTML Documentation with FORD

FORD (FORtran Documenter) generates beautiful, navigable HTML documentation from your Fortran source code.

**Step 1: Install FORD**
```bash
pip install ford
```

**Step 2: (Optional) Install Graphviz for Visual Diagrams**
```bash
# macOS
brew install graphviz

# Ubuntu/Debian
sudo apt-get install graphviz

# Windows
# Download from https://graphviz.org/download/
```

**Step 3: Generate Documentation**
```bash
# Navigate to project root
cd /path/to/lsdaks

# Generate documentation
ford ford.md
```

This creates a `doc/` directory with all HTML files.

**Step 4: View Documentation**
```bash
# macOS
open doc/index.html

# Linux
xdg-open doc/index.html

# Windows
start doc/index.html

# Or just open the file directly in your browser:
# file:///Users/guilherme.canella/Documents/lsdaks/doc/index.html
```

**What's Included:**
- ✅ **Module hierarchy** with interactive call graphs (if Graphviz installed)
- ✅ **Procedure documentation** with all parameters and return values
- ✅ **Source code browser** with syntax highlighting
- ✅ **Search functionality** to find functions/modules quickly
- ✅ **Mathematical formulas** rendered from LaTeX
- ✅ **Dependency diagrams** showing module relationships
- ✅ **Cross-references** between related code sections

**Tip:** Bookmark `doc/index.html` for quick access while coding!

### Project Documentation Files

- **README.md** (this file): User guide, installation, usage
- **CLAUDE.md**: Technical guidance for AI assistants
- **PROJECT_CONTEXT.md**: Detailed implementation notes (Portuguese)
- **MIXING_EQUIVALENCE.md**: Mixing convention between C++ and Fortran
- **ford.md**: FORD documentation generator configuration

## Physical Background

### The 1D Hubbard Model

The Hubbard Hamiltonian describes interacting electrons on a lattice:

$$
H = -t \sum_{\langle i,j \rangle, \sigma} (c^\dagger_{i\sigma} c_{j\sigma} + h.c.)
    + U \sum_i n_{i\uparrow} n_{i\downarrow}
    + \sum_{i\sigma} V^{\text{ext}}_i n_{i\sigma}
$$

- **Hopping term**: Kinetic energy (bandwidth ~ 4t)
- **Hubbard term**: On-site interaction (U > 0 repulsive, U < 0 attractive)
- **External potential**: Confining or disorder potentials

### Bethe Ansatz Solution

For the 1D case, the Bethe Ansatz provides exact eigenstates via the Lieb-Wu equations:

$$
e^{ik_j L} \prod_{\alpha=1}^M \frac{k_j - \Lambda_\alpha + iU/2}{k_j - \Lambda_\alpha - iU/2} = 1
$$

$$
\prod_{j=1}^N \frac{\Lambda_\alpha - k_j + iU/2}{\Lambda_\alpha - k_j - iU/2}
= \prod_{\beta \neq \alpha} \frac{\Lambda_\alpha - \Lambda_\beta + iU}{\Lambda_\alpha - \Lambda_\beta - iU}
$$

These are solved numerically using Newton-Raphson with analytical Jacobian.

### DFT-LSDA Framework

The many-body ground state energy is computed via density functional theory:

$$
E[n_\uparrow, n_\downarrow] = T_s[n] + E_{\text{Hartree}}[n] + E_{xc}[n_\uparrow, n_\downarrow] + \int V_{\text{ext}}(r) n(r) dr
$$

The exchange-correlation functional `E_xc` is obtained from Bethe Ansatz solutions and interpolated with bicubic splines in density and magnetization.

## Performance

### Bethe Ansatz Solver

The implemented solver is Newton-Raphson with an analytical Jacobian and continuation in U. Broyden is not implemented; timing depends on the grid and convergence history.

### SCF Convergence

Typical convergence in 50-200 iterations depending on:
- **Mixing parameter**: `α = 0.05` is conservative and stable
- **Adaptive mixing**: Helps with difficult cases
- **Potential type**: Smooth potentials converge faster

## Validation

### Energy Accuracy

Measured on 2026-09-20 using `build/cpp/lsdaks_cpp`, OBC, U=4 and the native U=4 table:

| Potential | C++ E/L | Fortran E/L | abs. difference | status |
|---|---:|---:|---:|---|
| uniform, L=90, 45/45, tolerances 1e-10 | -0.565718185 | -0.565718185262 | below the 5e-10 print resolution | reproduced; C++ output kept at `build/cpp/ref_uniform_u4` |
| harmonic `k=0.02`, L=20, 2/2 | -0.306692128 | -0.306692128394 | below the 5e-10 print resolution | historical; reproduced ad hoc on 2026-09-20 by rebuilding the input, no versioned input or script |
| double barrier `(3,3,-3,20)`, L=20, 2/2 | -0.971707225 | -0.971707224930 | below the 5e-10 print resolution | historical; reproduced ad hoc on 2026-09-20 by rebuilding the input, no versioned input or script |

The C++ prints the energy with `%12.9g` (`original/lsdaks.cc:44`), i.e. 9 significant digits,
so its last digit carries an uncertainty of ±5e-10. The raw differences
(`2.62e-10` for the uniform row, `3.94e-10` for the harmonic one, `7.0e-11` for
the double barrier) are all inside that print noise and must not be read as resolved
agreements at those magnitudes: 5e-10 is the floor of what this comparison can measure.
They are print resolution, not a measured level of agreement.

All three rows were measured on 2026-09-20 with the documented C++ executable
(`build/cpp/lsdaks_cpp`) and the potential mapping given in this README, and the harmonic and
double-barrier numbers agree with the values recorded earlier in the project. What separates
the rows is packaging, not physics: only the uniform case has its C++ output kept in the tree
(`build/cpp/ref_uniform_u4`). The last two rows have **no versioned input file and no
comparison script**, so reproducing them means rebuilding the inputs by hand from the
parameters in the table — which is exactly how they were checked. They are *historical* in
that narrow sense: reproducible in principle, just not from the repository alone. The legacy C++
type-4 impurity is not directly comparable: it
writes a six-site pattern and may exceed `1..L`; Fortran `impurity_single` means one physical
site.

### Known differences from C++

- The six Fortran potential-generator families are not a one-to-one equivalent of C++'s 13 types; Fortran has no `lsda_simetria.cc` path or interactive `r/m/s` loop.
- Intentionally not reproduced: out-of-range writes in C++ potential types 4 and 13, type-6's uninitialized final site, negative `Mix` from `DwMix`, non-refining `integral_1`, and `TOL=1e-16`.

### Physics Validation

- ⚠️ **U=0 (free fermions)**: Analytical reference used by the test suite only. The
  executable cannot run this case at all (see "An XC table is mandatory" above), so it
  is **not** validated end to end.
- ✅ **Half-filling (n=1)**: Matches Essler et al. reference values
- ✅ **Particle conservation**: `∫n dx = N` with error < 1e-12
- ⚠️ **Energy functional**: Stationary at the fixed point only up to the mismatch between
  `V_xc^σ` and `∂e_xc/∂n_σ` (~5e-4 away from `n = 1`, **no small bound at `n = 1`**),
  because `e_xc` and `V_xc` come from separate
  splines over separate tabulated columns; see the `xc_smoothing_width` note above. The
  opt-in smoothed potential (`w > 0`) is additionally non-variational by a finite amount.

#### Internal consistency checks (not external validation)

- ✅ **Bethe Ansatz Jacobian**: analytical Jacobian vs. finite differences of the **same**
  residual, agreement < 1e-10. This checks that the Jacobian matches the residual as
  coded; it says nothing about the residual being the correct Lieb-Wu equation. A missing
  `sin k` factor in the residual survived this check until phase 4 precisely because both
  sides used it. The residual is correct today, but the check itself is internal.

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Follow Fortran coding conventions (see CLAUDE.md)
4. Add tests for new features
5. Ensure all tests pass (`fpm test`)
6. Document code with FORD-style comments
7. Submit a pull request

### Coding Conventions

- **Modules**: `snake_case` (e.g., `bethe_equations`)
- **Types**: `snake_case_t` (e.g., `system_params_t`)
- **Functions**: `snake_case` (e.g., `solve_newton`)
- **Constants**: `UPPER_SNAKE_CASE` (e.g., `ITER_MAX`)
- **Precision**: Always use `real(dp)` from `lsda_constants`
- **Documentation**: FORD-compliant docstrings

## Citation

If you use this code in your research, please cite:

```bibtex
@software{lsda_hubbard_fortran,
  author = {Canella, Guilherme},
  title = {LSDA-Hubbard-Fortran: Local Spin Density Approximation for the 1D Hubbard Model},
  year = {2025},
  url = {https://github.com/yourusername/lsdaks},
  version = {0.1.0}
}
```

And the original Bethe Ansatz paper:

```bibtex
@article{lieb1968exact,
  title = {Absence of Mott Transition in an Exact Solution of the Short-Range, One-Band Model in One Dimension},
  author = {Lieb, Elliott H. and Wu, F. Y.},
  journal = {Phys. Rev. Lett.},
  volume = {20},
  pages = {1445--1448},
  year = {1968},
  doi = {10.1103/PhysRevLett.20.1445}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **E.H. Lieb and F.Y. Wu** for the Bethe Ansatz solution
- **Fortran Package Manager** developers for excellent build tooling
- **LAPACK/BLAS** community for high-performance linear algebra
- **Fortuno** developers for the modern Fortran testing framework

## Contact

**Guilherme Canella**
📧 guycanella@gmail.com
🐙 GitHub: [@gcanella](https://github.com/gcanella)

---

**Status**: tested with 20 registered suites / 1,455 source assertions | **License**: MIT 📄
