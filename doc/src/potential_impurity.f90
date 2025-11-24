!> Module for point impurity potentials
!!
!! Implements various types of point impurities:
!! - Single impurity at a fixed position
!! - Multiple impurities with individual amplitudes
!! - Random impurities with fixed concentration
!!
!! Physical interpretation:
!! - V_imp > 0: Repulsive impurity (barrier)
!! - V_imp < 0: Attractive impurity (potential well, can create bound states)
module potential_impurity
    use lsda_constants, only: dp
    use lsda_errors, only: ERROR_SUCCESS, ERROR_OUT_OF_BOUNDS, &
                           ERROR_SIZE_MISMATCH, ERROR_INVALID_CONCENTRATION
    implicit none
    private
    
    public :: potential_impurity_single
    public :: potential_impurity_multiple
    public :: potential_impurity_random

contains

    !> Single point impurity: V(i) = V_imp if i = i_imp, else 0
    !!
    !! Creates a point impurity at a specific lattice site.
    !!
    !! @param[in]  V_imp    Impurity amplitude
    !! @param[in]  i_imp    Impurity position (1-indexed, must be in [1, L])
    !! @param[in]  L        Number of lattice sites
    !! @param[out] V        Potential array V(i) for i = 1..L
    !! @param[out] ierr     Error flag (ERROR_SUCCESS or ERROR_OUT_OF_BOUNDS)
    subroutine potential_impurity_single(V_imp, i_imp, L, V, ierr)
        real(dp), intent(in) :: V_imp
        integer, intent(in) :: i_imp, L
        real(dp), dimension(L), intent(out) :: V
        integer, intent(out) :: ierr

        ! Initialize
        V = 0.0_dp
        ierr = ERROR_SUCCESS

        ! Bounds check
        if (i_imp < 1 .or. i_imp > L) then
            ierr = ERROR_OUT_OF_BOUNDS
            return
        end if

        ! Set impurity
        V(i_imp) = V_imp
    end subroutine potential_impurity_single

    !> Multiple point impurities with individual amplitudes
    !!
    !! Creates multiple impurities, each with its own amplitude.
    !! If two impurities are at the same position, their amplitudes are added.
    !!
    !! @param[in]  V_imp_array   Array of impurity amplitudes (size N_imp)
    !! @param[in]  imp_positions Array of impurity positions (size N_imp, 1-indexed)
    !! @param[in]  L             Number of lattice sites
    !! @param[out] V             Potential array V(i) for i = 1..L
    !! @param[out] ierr          Error flag (ERROR_SUCCESS, ERROR_SIZE_MISMATCH, or ERROR_OUT_OF_BOUNDS)
    subroutine potential_impurity_multiple(V_imp_array, imp_positions, L, V, ierr)
        real(dp), intent(in) :: V_imp_array(:)
        integer, intent(in) :: imp_positions(:), L
        real(dp), dimension(L), intent(out) :: V
        integer, intent(out) :: ierr
        integer :: k, N_imp, pos

        ! Initialize
        V = 0.0_dp
        ierr = ERROR_SUCCESS

        ! Check array sizes match
        N_imp = size(imp_positions)
        if (size(V_imp_array) /= N_imp) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        ! Add each impurity
        do k = 1, N_imp
            pos = imp_positions(k)
            
            ! Bounds check
            if (pos < 1 .or. pos > L) then
                ierr = ERROR_OUT_OF_BOUNDS
                return
            end if
            
            ! Add impurity (use += to handle overlapping impurities)
            V(pos) = V(pos) + V_imp_array(k)
        end do
    end subroutine potential_impurity_multiple

    !> Random impurities with fixed concentration
    !!
    !! Creates randomly distributed impurities with fixed amplitude.
    !! The number of impurities is determined by: N_imp = round(concentration * L / 100)
    !!
    !! Example: L = 100, concentration = 10 -> N_imp = 10 impurities (10%)
    !!
    !! @param[in]  V_imp         Impurity amplitude (same for all impurities)
    !! @param[in]  concentration Percentage of sites with impurities (0 < c <= 100)
    !! @param[in]  L             Number of lattice sites
    !! @param[in]  seed          Random seed for reproducibility (optional, use system time if < 0)
    !! @param[out] V             Potential array V(i) for i = 1..L
    !! @param[out] imp_positions Positions where impurities were placed (allocated inside)
    !! @param[out] ierr          Error flag (ERROR_SUCCESS or ERROR_INVALID_CONCENTRATION)
    subroutine potential_impurity_random(V_imp, concentration, L, seed, V, imp_positions, ierr)
        real(dp), intent(in) :: V_imp, concentration
        integer, intent(in) :: L, seed
        real(dp), dimension(L), intent(out) :: V
        integer, allocatable, intent(out) :: imp_positions(:)
        integer, intent(out) :: ierr
        integer :: N_imp, k, pos, attempt
        integer, allocatable :: seed_array(:)
        integer :: seed_size
        logical :: position_taken(L)
        real(dp) :: rand_val

        V = 0.0_dp
        ierr = ERROR_SUCCESS
        position_taken = .false.

        if (concentration <= 0.0_dp .or. concentration > 100.0_dp) then
            ierr = ERROR_INVALID_CONCENTRATION
            return
        end if

        N_imp = nint(concentration * real(L, dp) / 100.0_dp)
        if (N_imp == 0) N_imp = 1
        if (N_imp > L) N_imp = L 

        allocate(imp_positions(N_imp))

        if (seed >= 0) then
            call random_seed(size=seed_size)
            allocate(seed_array(seed_size))
            seed_array = seed
            call random_seed(put=seed_array)
            deallocate(seed_array)
        else
            call random_seed()
        end if

        ! Place impurities randomly (without replacement)
        do k = 1, N_imp
            attempt = 0
            do
                call random_number(rand_val)
                pos = 1 + int(rand_val * L)
                
                ! Make sure position is valid and not already taken
                if (pos >= 1 .and. pos <= L .and. .not. position_taken(pos)) then
                    position_taken(pos) = .true.
                    imp_positions(k) = pos
                    V(pos) = V_imp
                    exit
                end if
                
                ! Safety check: avoid infinite loop
                attempt = attempt + 1
                if (attempt > 1000 * L) then
                    ierr = ERROR_OUT_OF_BOUNDS 
                    return
                end if
            end do
        end do

        call sort_array(imp_positions)
    end subroutine potential_impurity_random

    !> Simple insertion sort for small arrays
    subroutine sort_array(arr)
        integer, intent(inout) :: arr(:)
        integer :: i, j, temp, n
        
        n = size(arr)
        do i = 2, n
            temp = arr(i)
            j = i - 1
            do while (j >= 1)
                if (arr(j) <= temp) exit
                arr(j + 1) = arr(j)
                j = j - 1
            end do
            arr(j + 1) = temp
        end do
    end subroutine sort_array

end module potential_impurity