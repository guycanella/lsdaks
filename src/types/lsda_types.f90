!> Derived types for LSDA-Hubbard calculations
!!
!! This module defines the core data structures used throughout the LSDA code.
!! These types encapsulate system parameters and state variables for DFT-LSDA
!! calculations on the 1D Hubbard model.
module lsda_types
    use lsda_constants, only: dp
    implicit none

    public :: system_params_t

    !> System parameters for LSDA-Hubbard calculations
    !!
    !! This type stores all the physical and numerical parameters that define
    !! a calculation for the 1D Hubbard model with an external potential.
    !!
    !! **Physical parameters:**
    !! - Number of sites, electrons (spin-up/down)
    !! - Hubbard interaction strength U
    !! - External potential type
    !!
    !! **Numerical parameters:**
    !! - Boundary conditions (open, periodic, twisted)
    !! - Convergence tolerance
    !! - Symmetry flags
    !!
    !! @note All parameters must be set before calling the KS cycle.
    !! @note For twisted BC (bc=2), the phase parameter must be provided.
    !!    
    type :: system_params_t
        integer :: Na           !< Number of lattice sites
        integer :: Nup          !< Number of spin-up electrons (N↑)
        integer :: Ndown        !< Number of spin-down electrons (N↓)
        integer :: bc           !< Boundary conditions: 0=open, 1=periodic, 2=twisted
        real(dp) :: U           !< Hubbard interaction strength (in units of hopping t=1)
        real(dp) :: phase       !< Twist angle for twisted BC (in units of π), only used if bc=2
    end type system_params_t
end module lsda_types