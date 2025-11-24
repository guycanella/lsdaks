module mixing_schemes
    use lsda_constants, only: dp
    use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT, &
                                            ERROR_SIZE_MISMATCH
    implicit none
    private

    public :: linear_mixing
    ! public :: anderson_mixing
    ! public :: broyden_mixing

contains
    !> @brief Linear density mixing for SCF convergence
    !!
    !! Implements: n_mixed = (1-α)·n_old + α·n_new
    !!
    !! Linear mixing damps oscillations in SCF cycles by combining
    !! input and output densities.
    !!
    !! Convention: α represents the weight of NEW density
    !! - α small (e.g., 0.05): conservative, keeps 95% old (like C++ Mix=0.95)
    !! - α large (e.g., 0.5): aggressive, 50/50 mix
    !!
    !! Default: α = 0.05 (equivalent to original C++ code with Mix=0.95)
    !!
    !! @param[in] n_new Output density from diagonalization (length L)
    !! @param[in] n_old Input density used to build Hamiltonian (length L)
    !! @param[in] alpha Mixing parameter (0 < α ≤ 1), weight of new density
    !! @param[out] n_mixed Mixed density for next iteration (length L)
    !! @param[in] L System size
    !! @param[out] ierr Error code (0 = success)
    subroutine linear_mixing(n_new, n_old, alpha, n_mixed, L, ierr)
        integer, intent(in) :: L
        real(dp), intent(in) :: n_new(:), n_old(:)
        real(dp), intent(in) :: alpha
        real(dp), intent(out) :: n_mixed(:)
        integer, intent(out) :: ierr

        integer :: i

        if (L <= 0) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (size(n_new) /= L .or. size(n_old) /= L .or. size(n_mixed) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        if (alpha <= 0.0_dp .or. alpha > 1.0_dp) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        ierr = ERROR_SUCCESS

        n_mixed = (1.0_dp - alpha) * n_old + alpha * n_new
    end subroutine linear_mixing
end module mixing_schemes