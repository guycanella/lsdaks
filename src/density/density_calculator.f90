module density_calculator
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use lsda_constants, only: dp
    use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT, &
                           ERROR_SIZE_MISMATCH, ERROR_UNPHYSICAL_DENSITY
    implicit none
    private

    real(dp), parameter :: TOL = 1.0e-10_dp

    public :: compute_occupations
    public :: compute_density_spin
    public :: compute_total_density
    public :: verify_particle_number
    public :: check_density_bounds

    interface compute_density_spin
        module procedure compute_density_spin_real
        module procedure compute_density_spin_complex
        module procedure compute_density_spin_occ_real
        module procedure compute_density_spin_occ_complex
    end interface compute_density_spin

contains
    !> @brief Occupation numbers of the Kohn-Sham levels at T = 0
    !!
    !! Fills the `n_elec` lowest levels, but shares the occupation of an OPEN
    !! (partially filled) degenerate shell equally among its members, instead
    !! of picking the `n_elec` first eigenvectors with weight 1.
    !!
    !! Why this matters: when level `n_elec` is degenerate with level
    !! `n_elec + 1` (PBC with a ±k pair and even N_σ, antiperiodic BC at θ = π,
    !! accidental degeneracies), integer occupation makes the density depend on
    !! whichever arbitrary basis LAPACK happened to return inside the degenerate
    !! subspace. With PBC and a uniform potential LAPACK returns the real
    !! cos/sin combinations of the k = ±m pair, and occupying only one of them
    !! creates a spurious period-2 density wave which feeds back into V_eff and
    !! turns the SCF cycle into a limit cycle that never converges.
    !!
    !! The shell is located exactly as the C++ reference does it
    !! (`update_degen`, original/lsdaks.cc:85-137): starting at `n_elec`, the
    !! window is widened while CONSECUTIVE NEIGHBOURS are closer than `deg_tol`.
    !! This is deliberately *not* the (simpler) criterion
    !! |eigvals(j) - eigvals(n_elec)| < deg_tol: the two agree on exact
    !! degeneracies but differ when the spectrum has a cascade of nearby levels,
    !! and the C++ behaviour is the one we must reproduce.
    !!
    !!   occ(j) = 1                                   for j < level_min
    !!   occ(j) = (n_elec - (level_min - 1)) / g      for level_min <= j <= level_max
    !!   occ(j) = 0                                   for j > level_max
    !!
    !! with g = level_max - level_min + 1 the shell degeneracy. Note the partial
    !! weight is UNIFORM over the whole block, and that sum(occ) = n_elec holds
    !! algebraically: no renormalisation of the density is needed afterwards.
    !!
    !! Unlike the C++ original, the window search is clamped to the available
    !! spectrum: `level_min - 1 >= 1` and `level_max + 1 <= size(eigvals)`. The
    !! C++ reads erg[0] / erg[num_eigen+1] out of bounds and relies on the
    !! `num_eigen = min(Ne + 5, Na)` slack to make it harmless.
    !!
    !! `n_elec == 0` (a fully polarised channel) is handled explicitly: all
    !! occupations are zero and `eigvals` is never indexed.
    !!
    !! @param[in]  eigvals Eigenvalues in ASCENDING order (length >= n_elec)
    !! @param[in]  n_elec  Number of electrons in this spin channel (>= 0)
    !! @param[in]  deg_tol Degeneracy line width (use DEG_TOL = 1.0e-10)
    !! @param[out] occ     Occupation of each level, same length as eigvals
    !! @param[out] ierr    Error code (0 = success)
    subroutine compute_occupations(eigvals, n_elec, deg_tol, occ, ierr)
        real(dp), intent(in) :: eigvals(:)
        integer, intent(in) :: n_elec
        real(dp), intent(in) :: deg_tol
        real(dp), intent(out) :: occ(:)
        integer, intent(out) :: ierr

        integer :: n_levels, level_min, level_max
        real(dp) :: g, partial_fill

        n_levels = size(eigvals)

        if (n_levels <= 0 .or. n_elec < 0 .or. n_elec > n_levels) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (.not. ieee_is_finite(deg_tol) .or. deg_tol < 0.0_dp) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (size(occ) /= n_levels) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        ierr = ERROR_SUCCESS
        occ(:) = 0.0_dp

        ! Fully polarised / empty channel: nothing to occupy, and eigvals(0)
        ! must not be touched (the C++ does touch it; see the note above).
        if (n_elec == 0) return

        ! Widen the Fermi shell downwards while consecutive neighbours coincide.
        level_min = n_elec
        do while (level_min > 1)
            if (abs(eigvals(level_min) - eigvals(level_min - 1)) >= deg_tol) exit
            level_min = level_min - 1
        end do

        ! ... and upwards.
        level_max = n_elec
        do while (level_max < n_levels)
            if (abs(eigvals(level_max + 1) - eigvals(level_max)) >= deg_tol) exit
            level_max = level_max + 1
        end do

        g = real(level_max - level_min + 1, dp)
        partial_fill = real(n_elec - (level_min - 1), dp) / g

        if (level_min > 1) occ(1:level_min - 1) = 1.0_dp
        occ(level_min:level_max) = partial_fill
    end subroutine compute_occupations

    !> @brief Compute density from real eigenvectors, integer occupation
    !!
    !! Calculates n_σ(i) = Σ_{j<=n_elec} |ψ_j(i)|² for occupied states at T=0.
    !!
    !! This variant fills the n_elec lowest columns with weight 1 and is only
    !! correct for a CLOSED Fermi shell. Production code (the SCF cycle) uses
    !! the occupation-array variant fed by `compute_occupations`, which also
    !! handles an open degenerate shell.
    !!
    !! `n_elec = 0` is legal and returns a zero density (fully polarised
    !! channel); only `n_elec < 0` is rejected.
    !!
    !! @param[in] eigvecs Real eigenvector matrix (L × L), columns are eigenvectors
    !! @param[in] L System size
    !! @param[in] n_elec Number of occupied electrons (for this spin), >= 0
    !! @param[out] density Electron density at each site (length L)
    !! @param[out] ierr Error code (0 = success)
    subroutine compute_density_spin_real(eigvecs, L, n_elec, density, ierr)
        real(dp), intent(in) :: eigvecs(:, :)
        integer, intent(in) :: L, n_elec
        real(dp), intent(out) :: density(:)
        integer, intent(out) :: ierr

        integer :: i

        if (L <= 0 .or. n_elec < 0 .or. n_elec > L) then
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

    !> @brief Compute density from complex eigenvectors, integer occupation
    !!
    !! Calculates n_σ(i) = Σ_{j<=n_elec} |ψ_j(i)|² for occupied states at T=0.
    !! Same caveat as the real version: correct only for a closed Fermi shell.
    !!
    !! @param[in] eigvecs Complex eigenvector matrix (L × L), columns are eigenvectors
    !! @param[in] L System size
    !! @param[in] n_elec Number of occupied electrons (for this spin), >= 0
    !! @param[out] density Real electron density at each site (length L)
    !! @param[out] ierr Error code (0 = success)
    subroutine compute_density_spin_complex(eigvecs, L, n_elec, density, ierr)
        complex(dp), intent(in) :: eigvecs(:, :)
        integer, intent(in) :: L, n_elec
        real(dp), intent(out) :: density(:)
        integer, intent(out) :: ierr

        integer :: i

        if (L <= 0 .or. n_elec < 0 .or. n_elec > L) then
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

    !> @brief Compute density from real eigenvectors with fractional occupation
    !!
    !! n_σ(i) = Σ_j occ(j) |ψ_j(i)|², summed only up to the last non-zero
    !! occupation. `occ` comes from `compute_occupations`, so an open degenerate
    !! Fermi shell contributes with a uniform fractional weight and the density
    !! no longer depends on the arbitrary basis LAPACK chose inside the shell.
    !!
    !! An all-zero `occ` (fully polarised channel) returns a zero density.
    !!
    !! @param[in] eigvecs Real eigenvector matrix (L × n_levels), columns are eigenvectors
    !! @param[in] L System size
    !! @param[in] occ Occupation of each level (length <= size(eigvecs, 2))
    !! @param[out] density Electron density at each site (length L)
    !! @param[out] ierr Error code (0 = success)
    subroutine compute_density_spin_occ_real(eigvecs, L, occ, density, ierr)
        real(dp), intent(in) :: eigvecs(:, :)
        integer, intent(in) :: L
        real(dp), intent(in) :: occ(:)
        real(dp), intent(out) :: density(:)
        integer, intent(out) :: ierr

        integer :: i, level_max

        if (L <= 0) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (size(eigvecs, 1) /= L .or. size(eigvecs, 2) < size(occ)) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        if (size(density) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        ierr = ERROR_SUCCESS

        level_max = last_occupied(occ)
        density(:) = 0.0_dp
        if (level_max < 1) return

        do i = 1, L
            density(i) = sum(occ(1:level_max) * eigvecs(i, 1:level_max)**2)
        end do
    end subroutine compute_density_spin_occ_real

    !> @brief Compute density from complex eigenvectors with fractional occupation
    !!
    !! Complex (twisted BC) counterpart of `compute_density_spin_occ_real`.
    !! Mirrors `update_dens_twist` in the C++ reference
    !! (original/lsda_twist.cc:13-54), which uses the same open-shell window and
    !! the same uniform partial weight.
    !!
    !! @param[in] eigvecs Complex eigenvector matrix (L × n_levels)
    !! @param[in] L System size
    !! @param[in] occ Occupation of each level (length <= size(eigvecs, 2))
    !! @param[out] density Real electron density at each site (length L)
    !! @param[out] ierr Error code (0 = success)
    subroutine compute_density_spin_occ_complex(eigvecs, L, occ, density, ierr)
        complex(dp), intent(in) :: eigvecs(:, :)
        integer, intent(in) :: L
        real(dp), intent(in) :: occ(:)
        real(dp), intent(out) :: density(:)
        integer, intent(out) :: ierr

        integer :: i, level_max

        if (L <= 0) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (size(eigvecs, 1) /= L .or. size(eigvecs, 2) < size(occ)) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        if (size(density) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        ierr = ERROR_SUCCESS

        level_max = last_occupied(occ)
        density(:) = 0.0_dp
        if (level_max < 1) return

        do i = 1, L
            density(i) = sum(occ(1:level_max) * &
                real(eigvecs(i, 1:level_max) * conjg(eigvecs(i, 1:level_max)), kind=dp))
        end do
    end subroutine compute_density_spin_occ_complex

    !> @brief Index of the highest level with non-zero occupation
    !!
    !! Returns 0 when every occupation is (numerically) zero. Occupations are
    !! compared against TOL rather than to 0 exactly, because a partial fill is
    !! a quotient of integers and an empty level may carry roundoff noise.
    !!
    !! @param[in] occ Occupation of each level
    !! @return Largest j with occ(j) > TOL, or 0 if there is none
    pure function last_occupied(occ) result(level_max)
        real(dp), intent(in) :: occ(:)
        integer :: level_max

        integer :: j

        level_max = 0
        do j = size(occ), 1, -1
            if (abs(occ(j)) > TOL) then
                level_max = j
                exit
            end if
        end do
    end function last_occupied

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