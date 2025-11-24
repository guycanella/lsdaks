module convergence_monitor
    use lsda_constants, only: dp
    use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT, &
                                            ERROR_SIZE_MISMATCH
    implicit none
    private

    type :: convergence_history_t
        integer :: max_iter
        integer :: current_iter
        real(dp), allocatable :: density_norms(:)
        real(dp), allocatable :: energies(:)
    end type convergence_history_t

    enum, bind(c)
        enumerator :: L1 = 1
        enumerator :: L2 = 2
        enumerator :: Linf = 3
    end enum

    public :: L1, L2, Linf
    public :: compute_density_difference
    public :: compute_density_norm
    public :: check_scf_convergence
    public :: update_convergence_history
    public :: convergence_history_t
    public :: init_convergence_history
    public :: cleanup_convergence_history

contains
    !> @brief Compute difference between new and old densities
    !!
    !! Calculates Δn(i) = n_new(i) - n_old(i)
    !!
    !! @param[in] n_new New density (length L)
    !! @param[in] n_old Old density (length L)
    !! @param[in] L System size
    !! @param[out] diff Density difference (length L)
    !! @param[out] ierr Error code (0 = success)
    subroutine compute_density_difference(n_new, n_old, L, diff, ierr)
        real(dp), intent(in) :: n_new(:), n_old(:)
        integer, intent(in) :: L
        real(dp), intent(out) :: diff(:)
        integer, intent(out) :: ierr

        if (L <= 0) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (size(n_new) /= L .or. size(n_old) /= L .or. size(diff) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        ierr = ERROR_SUCCESS
        diff = n_new - n_old
    end subroutine compute_density_difference

    !> @brief Compute norm of density difference
    !!
    !! Supports L1, L2, and L∞ norms:
    !! - L1:   ||Δn||₁ = Σᵢ |Δn(i)|
    !! - L2:   ||Δn||₂ = √(Σᵢ |Δn(i)|²)
    !! - L∞:   ||Δn||∞ = maxᵢ |Δn(i)|
    !!
    !! @param[in] delta_n Density difference (length L)
    !! @param[in] L System size
    !! @param[in] norm_type Type of norm (L1, L2, or Linf)
    !! @param[out] norm_value Computed norm
    !! @param[out] ierr Error code (0 = success)
    subroutine compute_density_norm(delta_n, L, norm_type, norm_value, ierr)
        real(dp), intent(in) :: delta_n(:)
        integer, intent(in) :: L
        integer, intent(in) :: norm_type
        real(dp), intent(out) :: norm_value
        integer, intent(out) :: ierr

        if (L <= 0) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (size(delta_n) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        ierr = ERROR_SUCCESS

        select case (norm_type)
        case (L1)
            norm_value = sum(abs(delta_n))
        case (L2)
            norm_value = sqrt(sum(delta_n**2))
        case (Linf)
            norm_value = maxval(abs(delta_n))
        case default
            ierr = ERROR_INVALID_INPUT
            norm_value = -1.0_dp
        end select
    end subroutine compute_density_norm

    !> @brief Check if SCF iteration has converged
    !!
    !! Convergence criterion: ||Δn||₂ < tol
    !!
    !! @param[in] delta_n Density difference (length L)
    !! @param[in] L System size
    !! @param[in] tol Convergence tolerance (default: 1.0e-6)
    !! @param[out] is_converged True if converged
    !! @param[out] ierr Error code (0 = success)
    subroutine check_scf_convergence(delta_n, L, tol, is_converged, ierr)
        real(dp), intent(in) :: delta_n(:)
        integer, intent(in) :: L
        real(dp), intent(in), optional :: tol
        logical, intent(out) :: is_converged
        integer, intent(out) :: ierr

        real(dp) :: norm_value
        real(dp) :: tolerance
        integer :: norm_type

        norm_type = L2

        if (present(tol)) then
            tolerance = tol
        else
            tolerance = 1.0e-6_dp
        end if

        call compute_density_norm(delta_n, L, norm_type, norm_value, ierr)

        if (ierr /= ERROR_SUCCESS) then
            is_converged = .false.
            return
        end if

        ierr = ERROR_SUCCESS
        is_converged = (norm_value < tolerance)
    end subroutine check_scf_convergence

    !> @brief Update convergence history with current iteration data
    !!
    !! @param[in] iteration Current iteration number
    !! @param[in] norm Density norm at this iteration
    !! @param[in] energy Total energy at this iteration
    !! @param[inout] history Convergence history object
    !! @param[out] ierr Error code (0 = success)
    subroutine update_convergence_history(iteration, norm, energy, history, ierr)
        integer, intent(in) :: iteration
        real(dp), intent(in) :: norm
        real(dp), intent(in) :: energy
        type(convergence_history_t), intent(inout) :: history
        integer, intent(out) :: ierr

        if (iteration < 1 .or. iteration > history%max_iter) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (.not. allocated(history%density_norms) .or. &
            .not. allocated(history%energies)) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        history%density_norms(iteration) = norm
        history%energies(iteration) = energy
        history%current_iter = iteration

        ierr = ERROR_SUCCESS
    end subroutine update_convergence_history

    !> @brief Initialize convergence history storage
    !!
    !! @param[out] history Convergence history object
    !! @param[in] max_iter Maximum number of iterations
    !! @param[out] ierr Error code (0 = success)
    subroutine init_convergence_history(history, max_iter, ierr)
        type(convergence_history_t), intent(out) :: history
        integer, intent(in) :: max_iter
        integer, intent(out) :: ierr
        
        if (max_iter <= 0) then
            ierr = ERROR_INVALID_INPUT
            return
        end if
        
        history%max_iter = max_iter
        history%current_iter = 0
        
        allocate(history%density_norms(max_iter))
        allocate(history%energies(max_iter))
        
        history%density_norms = 0.0_dp
        history%energies = 0.0_dp
        
        ierr = ERROR_SUCCESS
    end subroutine init_convergence_history

    !> @brief Deallocate convergence history
    !!
    !! @param[inout] history Convergence history object
    !! @param[out] ierr Error code (0 = success)
    subroutine cleanup_convergence_history(history, ierr)
        type(convergence_history_t), intent(inout) :: history
        integer, intent(out) :: ierr
        
        if (allocated(history%density_norms)) deallocate(history%density_norms)
        if (allocated(history%energies)) deallocate(history%energies)
        
        history%max_iter = 0
        history%current_iter = 0
        
        ierr = ERROR_SUCCESS
    end subroutine cleanup_convergence_history
end module convergence_monitor