!> Module for quasiperiodic potential (Aubry-André-Harper model)
!!
!! Implements the Aubry-André-Harper (AAH) quasiperiodic potential,
!! a paradigmatic model for studying localization transitions in 1D systems.
!!
!! Physical properties:
!! - Extended (delocalized) phase for λ < 2
!! - Critical phase at λ = 2 (multifractal wavefunctions)
!! - Localized phase for λ > 2 (exponentially localized states)
!!
!! The potential is given by:
!! V(i) = λ cos(2πβi + φ)
!!
!! where:
!! - λ (lambda) controls the potential strength
!! - β (beta) is typically an irrational number (golden ratio for maximum incommensurability)
!! - φ (phi) is a phase offset in radians
!!
!! Physics relevance:
!! - Models quasicrystalline structures
!! - Studies Anderson-like localization without disorder
!! - Relevant for cold atoms in bichromatic optical lattices
!! - Exhibits metal-insulator transitions
module potential_quasiperiodic
    use lsda_constants, only: dp, PI, TWOPI
    implicit none
    private

    public :: apply_potential_quasiperiodic

    !> Golden ratio (default for maximum incommensurability)
    real(dp), parameter :: GOLDEN_RATIO = 0.5_dp * (sqrt(5.0_dp) - 1.0_dp)

contains

    !> Quasiperiodic potential: V(i) = λ cos(2πβi + φ)
    !!
    !! Computes the Aubry-André-Harper potential at all lattice sites.
    !! For β = golden ratio, the system exhibits:
    !! - λ < 2: Extended phase (delocalized states)
    !! - λ = 2: Critical phase (multifractal wavefunctions)
    !! - λ > 2: Localized phase (exponentially localized states)
    !!
    !! @param[in]  lambda  Potential strength (typical range: 0 to 10)
    !! @param[in]  beta    Quasiperiodic frequency (default: golden ratio ≈ 0.618)
    !! @param[in]  phi     Phase offset in radians (default: 0)
    !! @param[in]  L       Number of lattice sites
    !! @param[out] V       Potential array V(i) for i = 1..L
    !! @param[out] ierr    Error flag
    !!
    !! @note Site indexing starts at 1 (Fortran convention).
    !!       For site i, the argument is 2πβ(i-1) + φ to start at i=0 physically.
    !! @note For β = golden ratio:
    !!       - λ < 2: All states are extended (delocalized)
    !!       - λ = 2: Critical phase (multifractal wavefunctions)
    !!       - λ > 2: All states are localized (exponentially decaying)
    subroutine apply_potential_quasiperiodic(lambda, beta, phi, L, V, ierr)
        use lsda_errors, only: ERROR_SUCCESS, ERROR_NEGATIVE_VALUE, ERROR_OUT_OF_BOUNDS
        real(dp), intent(in) :: lambda
        real(dp), intent(in) :: beta
        real(dp), intent(in) :: phi
        integer, intent(in) :: L
        real(dp), dimension(L), intent(out) :: V
        integer, intent(out) :: ierr

        integer :: i
        real(dp) :: argument

        ! Initialize
        V = 0.0_dp
        ierr = ERROR_SUCCESS

        ! Validate lambda (potential strength)
        if (lambda < 0.0_dp) then
            ierr = ERROR_NEGATIVE_VALUE
            return
        end if

        ! Validate beta (quasiperiodic frequency)
        if (beta <= 0.0_dp) then
            ierr = ERROR_NEGATIVE_VALUE
            return
        end if

        ! Validate L
        if (L < 1) then
            ierr = ERROR_OUT_OF_BOUNDS
            return
        end if

        ! Calculate potential at each site
        ! Use i-1 to make the physical site index start at 0
        do i = 1, L
            argument = TWOPI * beta * real(i - 1, dp) + phi
            V(i) = lambda * cos(argument)
        end do

    end subroutine apply_potential_quasiperiodic

end module potential_quasiperiodic
