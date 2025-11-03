!> Derived types for LSDA-Hubbard calculations
module lsda_types
    use lsda_constants, only: dp
    implicit none

    type :: system_params_t
        integer :: Na
        integer :: Nup
        integer :: Ndown
        integer :: bc
        real(dp) :: U
        real(dp) :: phase
    end type system_params_t
end module lsda_types