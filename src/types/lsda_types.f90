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
        integer :: L            !< Number of lattice sites
        integer :: Nup          !< Number of spin-up electrons (N↑)
        integer :: Ndown        !< Number of spin-down electrons (N↓)
        integer :: bc           !< Boundary conditions: 0=open, 1=periodic, 2=twisted
        real(dp) :: U           !< Hubbard interaction strength (in units of hopping t=1)
        !> Twist angle for twisted BC, **in radians**, only used if bc=2.
        !!
        !! Must lie in [0, 2π): this is the unit and the range that
        !! `apply_boundary_conditions_complex` and `validate_bc_parameters`
        !! expect. The user-facing input (`input_params_t%phase`, the `phase`
        !! key of the `&system` namelist and the `--phase` flag) is given in
        !! units of π, following the C++ reference, and is multiplied by π
        !! exactly once, in `convert_to_system_params`, when this field is
        !! filled. Anything that fills this field by hand must therefore give
        !! radians (e.g. 0.5 means half a radian, not π/2).
        real(dp) :: phase
    end type system_params_t
end module lsda_types