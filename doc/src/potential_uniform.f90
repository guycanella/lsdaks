!> Module for uniform (constant) external potential
!!
!! Implements a constant potential V(i) = V₀ for all sites.
!! This is equivalent to a global energy shift and does not affect
!! the physics of the system (only absolute energy scale).
module potential_uniform
    use lsda_constants, only: dp
    implicit none
    private

    public :: apply_potential_uniform

contains

    !> Uniform potential: V(i) = V₀ for all i
    !!
    !! Creates a constant potential across all lattice sites.
    !! This corresponds to a global energy shift.
    !!
    !! @param[in]  V0    Potential amplitude (constant value)
    !! @param[in]  L     Number of lattice sites
    !! @param[out] V     Potential array V(i) = V0 for i = 1..L
    !! @param[out] ierr  Error flag (always ERROR_SUCCESS for this potential)
    subroutine apply_potential_uniform(V0, L, V, ierr)
        use lsda_errors, only: ERROR_SUCCESS
        real(dp), intent(in) :: V0
        integer, intent(in) :: L
        real(dp), dimension(L), intent(out) :: V
        integer, intent(out) :: ierr

        V = V0
        ierr = ERROR_SUCCESS
    end subroutine apply_potential_uniform
end module potential_uniform