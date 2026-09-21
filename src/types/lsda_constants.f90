!> Physical and numerical constants for LSDA-Hubbard calculations

module lsda_constants
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    integer, parameter :: dp = real64                     ! Double precision (64 bits, ~16 digits)
    integer, parameter :: sp = int32                      ! Single precision (32 bits, ~8 digits)

    real(dp), parameter :: PI = 4.0_dp * atan(1.0_dp)     ! π = 3.14159...
    real(dp), parameter :: TWOPI = 2.0_dp * PI            ! 2π = 6.28318...

    real(dp), parameter :: U_SMALL = 1.0e-9_dp            ! Check if U ≈ 0
    real(dp), parameter :: TOL_DEFAULT = 1.0e-16_dp       ! Convergence tol

    ! Half-width of the numerical snap onto the Mott boundary n = 1.  Both
    ! the table generator's point evaluator and the spline consumer use this
    ! value, so a density in this round-off band is assigned the same
    ! one-sided limit before and after table generation.
    real(dp), parameter :: HALF_FILLING_SNAP_TOL = 1.0e-12_dp

    integer, parameter :: ITER_MAX = 10000                ! Max iterations SCF

    ! Mixing parameters (convention from original C++ code):
    ! C++ uses: V_new = Mix * V_old + (1-Mix) * V_calc (Mix = 0.95, keeps 95% old)
    ! Fortran uses: n_new = alpha * n_calc + (1-alpha) * n_old (alpha = weight of new)
    ! Therefore: alpha = 1 - Mix_cpp to maintain equivalence
    real(dp), parameter :: INITIAL_MIX = 0.95_dp          ! C++ mixing factor (keeps 95% old)
    real(dp), parameter :: MIX_ALPHA = 1.0_dp - INITIAL_MIX  ! Fortran mixing (5% new, 95% old)

    ! Floor for the adaptive mixing weight. Below this value each SCF step moves
    ! the effective potential by a negligible amount, so the cycle stops making
    ! progress while ||delta_n|| keeps shrinking proportionally to alpha: that is
    ! the false-convergence mechanism the potential residual criterion replaces.
    real(dp), parameter :: MIX_ALPHA_MIN = 0.005_dp       ! Minimum admissible alpha

    ! Half width of the line used to decide that two consecutive Kohn-Sham
    ! eigenvalues belong to the same degenerate shell. Two neighbours closer
    ! than DEG_TOL are treated as one level and share the occupation of the
    ! open Fermi shell equally (see density_calculator::compute_occupations).
    ! Same value as LINEWIDTH_ in the C++ reference (original/lsdaks.h:5).
    real(dp), parameter :: DEG_TOL = 1.0e-10_dp           ! Degeneracy line width

    ! Upper edge of the CONTINUOUS degeneracy detection (T20). Two neighbours
    ! closer than DEG_TOL are fully degenerate (link weight 1); farther apart
    ! than DEG_TOL_UPPER they are fully split (weight 0); in between the weight
    ! goes smoothly from 1 to 0, so the occupation numbers, the density and
    ! the band energy are continuous functions of the spectrum. The C++ uses a
    ! hard step at LINEWIDTH_ = 1e-10 (see compute_occupations for why that is
    ! a deliberate divergence). This is an ABSOLUTE energy width in units of t,
    ! not a relative tolerance: in very large PBC systems (roughly L >= 1e4),
    ! physically distinct levels can fall inside this transition interval.
    real(dp), parameter :: DEG_TOL_UPPER = 1.0e-6_dp      ! Degeneracy transition upper edge

    real(dp), parameter :: NEWTON_TOL = 1.0e-10_dp        ! Newton convergence tol
    integer, parameter :: NEWTON_MAX_ITER = 50            ! Max Newton iterations

    real(dp), parameter :: SCF_ENERGY_TOL = 1.0e-8_dp     ! SCF energy convergence tol
    real(dp), parameter :: SCF_DENSITY_TOL = 1.0e-6_dp    ! SCF density convergence tol
    real(dp), parameter :: SCF_POTENTIAL_TOL = 1.0e-6_dp  ! SCF potential residual convergence tol
end module lsda_constants
