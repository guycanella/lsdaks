module density_calculator
    use lsda_constants, only: dp
    use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT, &
                           ERROR_SIZE_MISMATCH, ERROR_UNPHYSICAL_DENSITY
    implicit none
    private

    real(dp), parameter :: TOL = 1.0e-10_dp

    public :: compute_density_spin
    public :: compute_total_density
    public :: verify_particle_number
    public :: check_density_bounds

    interface compute_density_spin
        module procedure compute_density_spin_real
        module procedure compute_density_spin_complex
    end interface compute_density_spin

contains
    !> @brief Compute density from real eigenvectors (Open/Periodic BC)
    !!
    !! Calculates n_σ(i) = Σ_j |ψ_j(i)|² for occupied states at T=0.
    !!
    !! @param[in] eigvecs Real eigenvector matrix (L × L), columns are eigenvectors
    !! @param[in] L System size
    !! @param[in] n_elec Number of occupied electrons (for this spin)
    !! @param[out] density Electron density at each site (length L)
    !! @param[out] ierr Error code (0 = success)
    subroutine compute_density_spin_real(eigvecs, L, n_elec, density, ierr)
        real(dp), intent(in) :: eigvecs(:, :)
        integer, intent(in) :: L, n_elec
        real(dp), intent(out) :: density(:)
        integer, intent(out) :: ierr

        integer :: i

        if (L <= 0 .or. n_elec <= 0 .or. n_elec > L) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (size(eigvecs, 1) /= L .or. size(eigvecs, 2) < n_elec) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        if (size(density) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        ierr = ERROR_SUCCESS

        do i = 1, L
            density(i) = sum(eigvecs(i, 1:n_elec)**2)
        end do
    end subroutine

    !> @brief Compute density from complex eigenvectors (Twisted BC)
    !!
    !! Calculates n_σ(i) = Σ_j |ψ_j(i)|² for occupied states at T=0.
    !!
    !! @param[in] eigvecs Complex eigenvector matrix (L × L), columns are eigenvectors
    !! @param[in] L System size
    !! @param[in] n_elec Number of occupied electrons (for this spin)
    !! @param[out] density Real electron density at each site (length L)
    !! @param[out] ierr Error code (0 = success)
    subroutine compute_density_spin_complex(eigvecs, L, n_elec, density, ierr)
        complex(dp), intent(in) :: eigvecs(:, :)
        integer, intent(in) :: L, n_elec
        real(dp), intent(out) :: density(:)
        integer, intent(out) :: ierr

        integer :: i

        if (L <= 0 .or. n_elec <= 0 .or. n_elec > L) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (size(eigvecs, 1) /= L .or. size(eigvecs, 2) < n_elec) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        if (size(density) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        ierr = ERROR_SUCCESS

        do i = 1, L
            density(i) = sum(real(eigvecs(i, 1:n_elec) * conjg(eigvecs(i, 1:n_elec)), kind=dp))
        end do
    end subroutine

    !> @brief Compute total density from spin-resolved densities
    !!
    !! Calculates n(i) = n_↑(i) + n_↓(i) at each site.
    !!
    !! @param[in] density_up Spin-up density (length L)
    !! @param[in] density_dw Spin-down density (length L)
    !! @param[out] density_total Total density n(i) (length L)
    !! @param[out] ierr Error code (0 = success)
    subroutine compute_total_density(density_up, density_dw, density_total, ierr)
        real(dp), intent(in) :: density_up(:), density_dw(:)
        real(dp), intent(out) :: density_total(:)
        integer, intent(out) :: ierr

        real(dp) :: n_total

        if (size(density_up) /= size(density_dw) .or. &
            size(density_up) /= size(density_total)) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        ierr = ERROR_SUCCESS

        density_total = density_up + density_dw
    end subroutine

    !> @brief Verify particle number conservation
    !!
    !! Checks if Σᵢ n_σ(i) = N_σ within tolerance TOL = 1e-10.
    !!
    !! @param[in] density Electron density (spin-up, spin-down, or total)
    !! @param[in] L System size
    !! @param[in] n_expected Expected number of particles
    !! @param[out] is_conserved True if |Σn - N| < TOL
    !! @param[out] ierr Error code (0 = success)
    subroutine verify_particle_number(density, L, n_expected, is_conserved, ierr)
        real(dp), intent(in) :: density(:)
        integer, intent(in) :: L, n_expected
        logical, intent(out) :: is_conserved
        integer, intent(out) :: ierr

        real(dp) :: n_computed

        if (L <= 0) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (size(density) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        ierr = ERROR_SUCCESS
        n_computed = sum(density)
        is_conserved = (abs(n_computed - real(n_expected, kind=dp)) < TOL)
    end subroutine

    !> @brief Check if densities satisfy physical bounds
    !!
    !! Verifies: 0 ≤ n_σ(i) and 0 ≤ n(i) = n_up(i) + n_dw(i) ≤ 2
    !!
    !! @param[in] density_up Spin-up density (length L)
    !! @param[in] density_dw Spin-down density (length L)
    !! @param[in] L System size
    !! @param[out] all_valid True if all densities are physical
    !! @param[out] ierr Error code (0 = success)
    subroutine check_density_bounds(density_up, density_dw, L, all_valid, ierr)
        real(dp), intent(in) :: density_up(:), density_dw(:)
        integer, intent(in) :: L
        logical, intent(out) :: all_valid
        integer, intent(out) :: ierr

        real(dp) :: n_total
        integer :: j

        if (L <= 0) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (size(density_up) /= L .or. size(density_dw) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        ierr = ERROR_SUCCESS
        all_valid = .true.

        do j = 1, L
            if (density_up(j) < -TOL .or. density_dw(j) < -TOL) then
                all_valid = .false.
                ierr = ERROR_UNPHYSICAL_DENSITY
                return
            end if

            n_total = density_up(j) + density_dw(j)
            if (n_total > 2.0_dp + TOL) then
                all_valid = .false.
                ierr = ERROR_UNPHYSICAL_DENSITY
                return
            end if
        end do
    end subroutine
end module density_calculator