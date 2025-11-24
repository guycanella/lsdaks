!> Module for barrier potentials
!!
!! Implements rectangular barrier potentials:
!! - Single barrier: One rectangular barrier
!! - Double barrier: Two barriers creating a quantum well
!!
!! Physical interpretation:
!! - V_bar > 0: Potential barrier (particle reflection)
!! - V_bar < 0: Potential well (particle attraction)
!! - Models quantum tunneling, resonant transmission
!! - Double barrier: Fabry-Pérot resonances, quasi-bound states
!!
!! Reference: Quantum mechanics textbooks (Griffiths, Cohen-Tannoudji)
module potential_barrier
    use lsda_constants, only: dp
    use lsda_errors, only: ERROR_SUCCESS, ERROR_OUT_OF_BOUNDS
    implicit none
    private
    
    public :: potential_barrier_single
    public :: potential_barrier_double

contains

    !> Single rectangular barrier: V(i) = V_bar for i_start ≤ i ≤ i_end, else 0
    !!
    !! Creates a rectangular barrier in a specified region.
    !!
    !! Physical properties:
    !! - V_bar > 0: Barrier (reflection, tunneling through barrier)
    !! - V_bar < 0: Well (attraction, bound states possible)
    !! - Width w = i_end - i_start + 1 determines tunneling probability
    !! - For V_bar >> E: Exponential decay of wavefunction in barrier
    !!
    !! Applications:
    !! - Quantum tunneling experiments
    !! - Scanning tunneling microscopy (STM)
    !! - Tunnel diodes
    !!
    !! @param[in]  V_bar    Barrier height (can be positive or negative)
    !! @param[in]  i_start  Starting position of barrier (1-indexed)
    !! @param[in]  i_end    Ending position of barrier (1-indexed)
    !! @param[in]  L        Number of lattice sites
    !! @param[out] V        Potential array V(i) for i = 1..L
    !! @param[out] ierr     Error flag (ERROR_SUCCESS or ERROR_OUT_OF_BOUNDS)
    !!
    !! @note Width: w = i_end - i_start + 1
    !! @note For w → 0: No barrier (delta-like)
    !! @note For w → ∞: Infinite barrier (no tunneling)
    subroutine potential_barrier_single(V_bar, i_start, i_end, L, V, ierr)
        real(dp), intent(in) :: V_bar
        integer, intent(in) :: i_start, i_end, L
        real(dp), dimension(L), intent(out) :: V
        integer, intent(out) :: ierr

        V = 0.0_dp
        ierr = ERROR_SUCCESS

        if (i_start < 1 .or. i_start > L) then
            ierr = ERROR_OUT_OF_BOUNDS
            return
        end if
        
        if (i_end < 1 .or. i_end > L) then
            ierr = ERROR_OUT_OF_BOUNDS
            return
        end if

        if (i_end < i_start) then
            ierr = ERROR_OUT_OF_BOUNDS
            return
        end if

        V(i_start:i_end) = V_bar
    end subroutine potential_barrier_single

    !> Double rectangular barrier: Two barriers creating a quantum well
    !!
    !! Creates two barriers separated by a quantum well, exactly matching C++ implementation.
    !! This geometry exhibits Fabry-Pérot-like resonances.
    !!
    !! C++ signature: double_barrier(int Na, double Vb, double Lb, double Vwell, double Lwell, double **v_ext)
    !!
    !! Geometry (centered at x0):
    !!   x_2         x_1        x0        x1         x2
    !!    |----Lb----|---Lwell/2---|---Lwell/2---|----Lb----|
    !!   [  Barreira  ][ Poço (Vwell) ][  Barreira  ]
    !!
    !! Physical properties:
    !! - Creates quasi-bound states in the well between barriers
    !! - Resonant tunneling at specific energies
    !! - Transmission peaks when E matches well energy levels
    !! - Well region typically has V_well < 0 (attractive potential)
    !!
    !! Applications:
    !! - Resonant tunneling diodes (RTD)
    !! - Quantum cascade lasers
    !! - Electron wave interferometry
    !!
    !! @param[in]  V_bar     Barrier height (positive for barrier)
    !! @param[in]  L_bar     Barrier width (each barrier has this width)
    !! @param[in]  V_well    Well potential (typically negative for attractive well)
    !! @param[in]  L_well    Well width (distance between inner edges of barriers)
    !! @param[in]  L         Number of lattice sites
    !! @param[out] V         Potential array V(i) for i = 1..L
    !! @param[out] ierr      Error flag (ERROR_SUCCESS or ERROR_OUT_OF_BOUNDS)
    !!
    !! @note Matches C++ double_barrier implementation exactly (lsda_potential.cc lines 166-194)
    !! @note SMALL = 1.0e-10 tolerance for region boundaries
    !! @note For L_well small: Strong coupling, large splitting
    !! @note For L_well large: Weak coupling, sharp resonances
    subroutine potential_barrier_double(V_bar, L_bar, V_well, L_well, L, V, ierr)
        real(dp), intent(in) :: V_bar, L_bar, V_well, L_well
        integer, intent(in) :: L
        real(dp), dimension(L), intent(out) :: V
        integer, intent(out) :: ierr

        real(dp) :: x0, x1, x2, x_1, x_2
        integer :: i
        real(dp), parameter :: SMALL = 1.0e-10_dp

        ierr = ERROR_SUCCESS
        V = 0.0_dp

        ! Calculate center position (matches C++ logic)
        ! C++: x0 = (Na%2 ? (double(Na/2)) + 1.0 : (double(Na/2)) + 0.5)
        if (mod(L, 2) == 1) then
            ! Odd L: center at (L/2) + 1 = (L+1)/2
            x0 = real(L/2, dp) + 1.0_dp
        else
            ! Even L: center at (L/2) + 0.5
            x0 = real(L/2, dp) + 0.5_dp
        end if

        ! Calculate region boundaries (matches C++ lines 171-175)
        x1 = x0 + L_well / 2.0_dp        ! Right edge of well
        x2 = x1 + L_bar                   ! Right edge of right barrier

        x_1 = x0 - L_well / 2.0_dp       ! Left edge of well
        x_2 = x_1 - L_bar                 ! Left edge of left barrier

        ! Apply potential to each site (matches C++ lines 179-191)
        do i = 1, L
            if (x_2 - SMALL < real(i, dp) .and. real(i, dp) < x_1 + SMALL) then
                ! Left barrier region
                V(i) = V_bar
            else if (x1 - SMALL < real(i, dp) .and. real(i, dp) < x2 + SMALL) then
                ! Right barrier region
                V(i) = V_bar
            else if (x_1 + SMALL <= real(i, dp) .and. real(i, dp) <= x1 - SMALL) then
                ! Well region (between barriers)
                V(i) = V_well
            else
                ! Outside regions
                V(i) = 0.0_dp
            end if
        end do
    end subroutine potential_barrier_double

end module potential_barrier