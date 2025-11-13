!> Module for uniform (constant) external potential
!!
!! Implements a constant potential V(i) = V₀ for all sites.
!! This is equivalent to a global energy shift and does not affect
!! the physics of the system (only absolute energy scale).
module potential_uniform
    use lsda_constants, only: dp
    implicit none
    private

    public :: potential_uniform

contains

    !> Uniform potential: V(i) = V₀ for all i
    !!
    !! Creates a constant potential across all lattice sites.
    !! This corresponds to a global energy shift.
    !!
    !! @param[in]  V0    Potential amplitude (constant value)
    !! @param[in]  L     Number of lattice sites
    !! @return     V     Potential array V(i) = V0 for i = 1..L
    subroutine potential_uniform(V0, L, V)
        real(dp), intent(in) :: V0
        integer, intent(in) :: L
        real(dp), dimension(L), intent(out) :: V

        V = V0
    end subroutine potential_uniform
end module potential_uniform