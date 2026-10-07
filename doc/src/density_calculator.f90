module density_calculator
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use lsda_constants, only: dp, DEG_TOL_UPPER
    use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT, &
                           ERROR_SIZE_MISMATCH, ERROR_UNPHYSICAL_DENSITY
    implicit none
    private

    real(dp), parameter :: TOL = 1.0e-10_dp

    public :: compute_occupations
    public :: link_weight
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
    !> @brief Zero-temperature Aufbau occupations with continuous near-degeneracy smearing
    !!
    !! Fills the `n_elec` lowest levels, but shares the occupation of an OPEN
    !! (partially filled) degenerate shell equally among its members, instead
    !! of picking the `n_elec` first eigenvectors with weight 1.
    !!
    !! Why this matters: when level `n_elec` is degenerate with level
    !! `n_elec + 1` (PBC with a ±k pair and even N_σ, antiperiodic BC at θ = π,
    !! accidental degeneracies, tunnel-split doublets in a trap whose core is
    !! a band insulator), integer occupation makes the density depend on
    !! whichever arbitrary basis LAPACK happened to return inside the degenerate
    !! subspace. With PBC and a uniform potential LAPACK returns the real
    !! cos/sin combinations of the k = ±m pair, and occupying only one of them
    !! creates a spurious period-2 density wave which feeds back into V_eff and
    !! turns the SCF cycle into a limit cycle that never converges.
    !!
    !! **Shell detection.** The shell is located along the same CONSECUTIVE
    !! NEIGHBOUR chain as the C++ reference (`update_degen`,
    !! original/lsdaks.cc:85-137): starting at `n_elec`, the chain is extended
    !! link by link while neighbours are close. This is deliberately *not* the
    !! (simpler) criterion |eigvals(j) - eigvals(n_elec)| < deg_tol: the two
    !! agree on exact degeneracies but differ on a cascade of nearby levels.
    !!
    !! **Deliberate divergence from the C++ (T20).** The C++ decides each link
    !! with a HARD step, `|e(j+1) - e(j)| < LINEWIDTH_ = 1e-10`. That makes the
    !! Kohn-Sham map n[V_eff] DISCONTINUOUS: a gap fluctuating around 1e-10 by
    !! roundoff flips a doublet at the Fermi level between (1/2, 1/2) and
    !! (1, 0). When LAPACK returns the doublet in a LOCALISED basis (the two
    !! arms of a harmonic trap separated by a filled core), that flip moves a
    !! whole electron from one arm to the other, breaks the reflection symmetry
    !! n(i) = n(L+1-i), and the SCF cycle has no fixed point to converge to.
    !! Here each link carries a CONTINUOUS weight.  Thus this is an effective
    !! C¹ smearing over the finite energy interval [DEG_TOL, DEG_TOL_UPPER],
    !! rather than a strictly sharp T = 0 occupation at every finite gap.
    !!
    !!   s(Δ) = 1                       for Δ <  deg_tol
    !!   s(Δ) = (1 - x)^2 (1 + 2x)      for deg_tol <= Δ < deg_tol_hi,
    !!                                  x = (Δ - deg_tol)/(deg_tol_hi - deg_tol)
    !!   s(Δ) = 0                       for Δ >= deg_tol_hi
    !!
    !! (a C¹ smoothstep, zero slope at both edges), and level j is connected to
    !! the Fermi level with the product w_j of the link weights between them
    !! (w_{n_elec} = 1). The integer (Aufbau) filling f0_j = 1 for j <= n_elec,
    !! 0 otherwise, is then redistributed inside the shell:
    !!
    !!   occ(j) = (1 - w_j) f0_j + w_j P / W,   P = Σ_k w_k f0_k,  W = Σ_k w_k
    !!
    !! i.e. every level pools a fraction w_j of its Aufbau charge and takes back
    !! the share w_j / W of the pool. Properties (in exact arithmetic):
    !!   * Σ_j occ(j) = n_elec algebraically (no renormalisation needed);
    !!   * 0 <= occ(j) <= 1 (convex combination of f0_j and P/W, both in [0,1]);
    !!   * occ is a continuous function of the eigenvalues;
    !!   * when every w_j is 0 or 1 - all gaps either below deg_tol or above
    !!     deg_tol_hi, which is the case for every exact degeneracy - the
    !!     result is bit-for-bit IDENTICAL to the C++ block rule
    !!       occ(j) = 1 for j < level_min, (n_elec - level_min + 1)/g inside
    !!       the block of size g, 0 above,
    !!     so all C++-parity results are unchanged.
    !! The density built from these occupations is the trace of the projector
    !! on the shell times the shared weight when the gap closes, which is
    !! invariant under the basis LAPACK picks inside the shell: the localised
    !! doublet and the symmetric doublet give the same density. Fixing the
    !! detection (not the basis) is therefore what removes the discontinuity.
    !!
    !! Unlike the C++ original, the chain is clamped to the available
    !! spectrum: the C++ reads erg[0] / erg[num_eigen+1] out of bounds and
    !! relies on the `num_eigen = min(Ne + 5, Na)` slack to make it harmless.
    !!
    !! `n_elec == 0` (a fully polarised channel) is handled explicitly: all
    !! occupations are zero and `eigvals` is never indexed.
    !!
    !! @param[in]  eigvals       Eigenvalues in ASCENDING order (length >= n_elec)
    !! @param[in]  n_elec        Number of electrons in this spin channel (>= 0)
    !! @param[in]  deg_tol       Degeneracy line width below which two neighbours
    !!                           are fully degenerate (use DEG_TOL = 1.0e-10)
    !! @param[out] occ           Occupation of each level, same length as eigvals
    !! @param[out] ierr          Error code (0 = success)
    !! @param[in]  deg_tol_hi    Optional upper edge of the continuous transition
    !!                           (>= deg_tol). Default: max(DEG_TOL_UPPER, deg_tol).
    !!                           Passing deg_tol_hi = deg_tol recovers the hard
    !!                           C++ step exactly. (Named deg_tol_hi, not
    !!                           deg_tol_upper: Fortran is case-insensitive and a
    !!                           dummy called deg_tol_upper would shadow the
    !!                           constant DEG_TOL_UPPER.)
    !! @param[out] shell_open    Optional. True when the degenerate chain is
    !!                           still open at the last level of `eigvals`, i.e.
    !!                           the shell may extend past the spectrum supplied.
    !!                           Callers that pass a PARTIAL spectrum must grow
    !!                           the window and call again; with a full spectrum
    !!                           there is nothing above and the flag is moot.
    subroutine compute_occupations(eigvals, n_elec, deg_tol, occ, ierr, deg_tol_hi, shell_open)
        real(dp), intent(in) :: eigvals(:)
        integer, intent(in) :: n_elec
        real(dp), intent(in) :: deg_tol
        real(dp), intent(out) :: occ(:)
        integer, intent(out) :: ierr
        real(dp), intent(in), optional :: deg_tol_hi
        logical, intent(out), optional :: shell_open

        integer :: n_levels, j
        real(dp) :: tol_hi, pool, weight_sum, shared
        real(dp) :: w(size(eigvals))

        if (present(shell_open)) shell_open = .false.

        n_levels = size(eigvals)

        if (n_levels <= 0 .or. n_elec < 0 .or. n_elec > n_levels) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (.not. ieee_is_finite(deg_tol) .or. deg_tol < 0.0_dp) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (present(deg_tol_hi)) then
            if (.not. ieee_is_finite(deg_tol_hi) .or. deg_tol_hi < deg_tol) then
                ierr = ERROR_INVALID_INPUT
                return
            end if
            tol_hi = deg_tol_hi
        else
            tol_hi = max(DEG_TOL_UPPER, deg_tol)
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

        ! Connection weight of every level to the Fermi level: the product of
        ! the link weights along the chain of consecutive neighbours. The chain
        ! stops as soon as a link is fully open (weight 0); beyond it w = 0.
        w(:) = 0.0_dp
        w(n_elec) = 1.0_dp

        do j = n_elec - 1, 1, -1
            w(j) = w(j + 1) * link_weight(eigvals(j + 1) - eigvals(j), deg_tol, tol_hi)
            if (w(j) <= 0.0_dp) exit
        end do

        do j = n_elec + 1, n_levels
            w(j) = w(j - 1) * link_weight(eigvals(j) - eigvals(j - 1), deg_tol, tol_hi)
            if (w(j) <= 0.0_dp) exit
        end do

        ! Pool the Aufbau charge of the connected levels (f0 = 1 for j <= n_elec)
        ! and hand it back in proportion to the connection weights. W >= 1
        ! always, because w(n_elec) = 1.
        pool = sum(w(1:n_elec))
        weight_sum = sum(w(1:n_levels))
        shared = pool / weight_sum

        ! Written in convex-combination form so w = 1 returns `shared`
        ! bit-for-bit, matching the C++ block rule for exact degeneracies.
        occ(1:n_elec) = (1.0_dp - w(1:n_elec)) + w(1:n_elec) * shared
        if (n_elec < n_levels) occ(n_elec + 1:n_levels) = w(n_elec + 1:n_levels) * shared

        ! The chain is still open at the top of the supplied spectrum: there may
        ! be further members of this shell that were never diagonalized, in which
        ! case `weight_sum` is short and `shared` too large. Only the caller
        ! knows whether more levels exist, so report instead of deciding here;
        ! when `eigvals` is the full spectrum this flag is meaningless and must
        ! be ignored.
        if (present(shell_open)) shell_open = w(n_levels) > 0.0_dp
    end subroutine compute_occupations

    !> @brief Continuous "same shell" weight of two consecutive levels
    !!
    !! Returns 1 for a gap below `tol_lo`, 0 for a gap at or above `tol_hi`,
    !! and the C¹ smoothstep 1 - x²(3 - 2x), x = (Δ - tol_lo)/(tol_hi -
    !! tol_lo), in between. This expanded complement form is deliberately used
    !! instead of the algebraically equivalent `(1 - x)^2(1 + 2x)`: it fixes
    !! the IEEE rounding of the transition weight (see its regression test).
    !! With tol_hi <= tol_lo it degenerates to the hard
    !! step of the C++ (`|Δ| < LINEWIDTH_`).
    !!
    !! @param[in] delta  Gap between the two levels (sign irrelevant)
    !! @param[in] tol_lo Gap below which the levels are fully degenerate
    !! @param[in] tol_hi Gap at or above which the levels are fully split
    !! @return           Weight in [0, 1]
    pure function link_weight(delta, tol_lo, tol_hi) result(s)
        real(dp), intent(in) :: delta, tol_lo, tol_hi
        real(dp) :: s

        real(dp) :: gap, x

        gap = abs(delta)

        if (gap < tol_lo) then
            s = 1.0_dp
        else if (gap >= tol_hi) then
            s = 0.0_dp
        else
            x = (gap - tol_lo) / (tol_hi - tol_lo)
            s = 1.0_dp - x * x * (3.0_dp - 2.0_dp * x)
        end if
    end function link_weight

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
