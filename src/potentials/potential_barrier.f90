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
    !! Creates two barriers separated by a distance d, forming a quantum well.
    !! This geometry exhibits Fabry-Pérot-like resonances.
    !!
    !! Physical properties:
    !! - Creates quasi-bound states in the well between barriers
    !! - Resonant tunneling at specific energies
    !! - Transmission peaks when E matches well energy levels
    !! - Well width: d = i2_start - i1_end - 1
    !!
    !! Applications:
    !! - Resonant tunneling diodes (RTD)
    !! - Quantum cascade lasers
    !! - Electron wave interferometry
    !!
    !! @param[in]  V_bar     Barrier height (same for both barriers)
    !! @param[in]  i1_start  Starting position of first barrier
    !! @param[in]  i1_end    Ending position of first barrier
    !! @param[in]  i2_start  Starting position of second barrier
    !! @param[in]  i2_end    Ending position of second barrier
    !! @param[in]  L         Number of lattice sites
    !! @param[out] V         Potential array V(i) for i = 1..L
    !! @param[out] ierr      Error flag (ERROR_SUCCESS or ERROR_OUT_OF_BOUNDS)
    !!
    !! @note Well width: d = i2_start - i1_end - 1
    !! @note For d small: Strong coupling, large splitting
    !! @note For d large: Weak coupling, sharp resonances
    subroutine potential_barrier_double(V_bar, i1_start, i1_end, i2_start, i2_end, L, V, ierr)
        real(dp), intent(in) :: V_bar
        integer, intent(in) :: i1_start, i1_end, i2_start, i2_end, L
        real(dp), dimension(L), intent(out) :: V
        integer, intent(out) :: ierr

        V = 0.0_dp
        ierr = ERROR_SUCCESS

        ! Bounds check for first barrier
        if (i1_start < 1 .or. i1_start > L) then
            ierr = ERROR_OUT_OF_BOUNDS
            return
        end if
        
        if (i1_end < 1 .or. i1_end > L) then
            ierr = ERROR_OUT_OF_BOUNDS
            return
        end if

        if (i1_end < i1_start) then
            ierr = ERROR_OUT_OF_BOUNDS
            return
        end if

        ! Bounds check for second barrier
        if (i2_start < 1 .or. i2_start > L) then
            ierr = ERROR_OUT_OF_BOUNDS
            return
        end if
        
        if (i2_end < 1 .or. i2_end > L) then
            ierr = ERROR_OUT_OF_BOUNDS
            return
        end if

        if (i2_end < i2_start) then
            ierr = ERROR_OUT_OF_BOUNDS
            return
        end if

        ! Check that barriers don't overlap
        if (i2_start <= i1_end) then
            ierr = ERROR_OUT_OF_BOUNDS
            return
        end if

        ! Set both barriers
        V(i1_start:i1_end) = V_bar
        V(i2_start:i2_end) = V_bar
    end subroutine potential_barrier_double

end module potential_barrier