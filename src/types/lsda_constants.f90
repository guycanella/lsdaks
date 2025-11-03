!> Physical and numerical constants for LSDA-Hubbard calculations

module lsda_constants
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    integer, parameter :: dp = real64                   ! Double precision (64 bits, ~16 digits)
    integer, parameter :: sp = int32                    ! Single precision (32 bits, ~8 digits)

    real(dp), parameter :: PI = 4.0_dp * atan(1.0_dp)   ! π = 3.14159...
    real(dp), parameter :: TWOPI = 2.0_dp * PI          ! 2π = 6.28318...

    real(dp), parameter :: SMALL = 1.0e-9_dp            ! Check if U ≈ 0
    real(dp), parameter :: TOL_DEFAULT = 1.0e-16_dp     ! Convergence tol

    integer, parameter :: ITER_MAX = 10000              ! Max iterations SCF
    real(dp), parameter :: INITIAL_MIX = 0.95_dp        ! Initial mixing factor
end module lsda_constants
