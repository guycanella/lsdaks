!> Module for harmonic trap potential
!!
!! Implements a parabolic confining potential centered at the middle of the chain.
!! Models optical traps in cold atom systems and creates shell structure.
module potential_harmonic
    use lsda_constants, only: dp
    implicit none
    private

    public :: potential_harmonic

contains

    !> Harmonic trap potential: V(i) = 0.5 * k * (i - i_center)²
    !!
    !! Creates a parabolic confining potential centered at the middle of the chain.
    !! The center is at i_center = (L+1)/2, which works for both odd and even L.
    !!
    !! Physical properties:
    !! - Confines particles at the center
    !! - Creates shell structure (similar to Landau levels)
    !! - Density is maximum at center, decays at edges
    !! - Has parity symmetry: V(i) = V(L+1-i)
    !!
    !! @param[in]  k    Spring constant (trap strength), k > 0
    !! @param[in]  L    Number of lattice sites
    !! @return     V    Potential array V(i) for i = 1..L
    !!
    !! @note For k → 0, reduces to uniform potential V = 0
    !! @note For k >> 1, creates strong confinement (Landau regime)
    subroutine potential_harmonic(k, L, V)
        real(dp), intent(in) :: k
        integer, intent(in) :: L
        real(dp), dimension(L), intent(out) :: V
        real(dp) :: i_center
        integer :: i

        i_center = real(L + 1, dp) / 2.0_dp

        V = [(0.5_dp * k * (real(i, dp) - i_center)**2, i = 1, L)]
    end subroutine potential_harmonic
end module potential_harmonic