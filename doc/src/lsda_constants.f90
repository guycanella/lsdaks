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

    integer, parameter :: ITER_MAX = 10000                ! Max iterations SCF

    ! Mixing parameters (convention from original C++ code):
    ! C++ uses: V_new = Mix * V_old + (1-Mix) * V_calc (Mix = 0.95, keeps 95% old)
    ! Fortran uses: n_new = alpha * n_calc + (1-alpha) * n_old (alpha = weight of new)
    ! Therefore: alpha = 1 - Mix_cpp to maintain equivalence
    real(dp), parameter :: INITIAL_MIX = 0.95_dp          ! C++ mixing factor (keeps 95% old)
    real(dp), parameter :: MIX_ALPHA = 1.0_dp - INITIAL_MIX  ! Fortran mixing (5% new, 95% old)

    real(dp), parameter :: NEWTON_TOL = 1.0e-10_dp        ! Newton convergence tol
    integer, parameter :: NEWTON_MAX_ITER = 50            ! Max Newton iterations

    real(dp), parameter :: SCF_ENERGY_TOL = 1.0e-8_dp     ! SCF energy convergence tol
    real(dp), parameter :: SCF_DENSITY_TOL = 1.0e-6_dp    ! SCF density convergence tol
end module lsda_constants
