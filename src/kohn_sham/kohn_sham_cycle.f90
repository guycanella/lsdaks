module kohn_sham_cycle
    use lsda_constants, only: dp, SCF_DENSITY_TOL, SCF_ENERGY_TOL, SCF_POTENTIAL_TOL, &
                              ITER_MAX, MIX_ALPHA
    use lsda_types, only: system_params_t
    use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT, &
                           ERROR_CONVERGENCE_FAILED, ERROR_SIZE_MISMATCH
    use boundary_conditions, only: apply_boundary_conditions, apply_boundary_conditions_complex
    use hamiltonian_builder, only: build_hamiltonian, build_hamiltonian_complex
    use lapack_wrapper, only: diagonalize_symmetric_real, diagonalize_hermitian_complex
    use density_calculator, only: compute_density_spin
    use xc_lsda, only: xc_lsda_t, get_vxc, get_exc
    use convergence_monitor, only: compute_density_difference, &
                                    convergence_history_t, init_convergence_history, &
                                    update_convergence_history, cleanup_convergence_history, &
                                    L2, compute_density_norm
    use mixing_schemes, only: linear_mixing
    use adaptive_mixing, only: adaptive_mix_t, adaptive_mix_init, adaptive_mix_update, &
                                adaptive_mix_get_alpha, adaptive_mix_reset
    implicit none
    private

    !> @brief Parameters for SCF cycle
    type :: scf_params_t
        integer :: max_iter = ITER_MAX             ! Maximum SCF iterations
        real(dp) :: density_tol = SCF_DENSITY_TOL  ! Diagnostic tolerance for ||Δn|| (no longer a convergence criterion)
        real(dp) :: energy_tol = SCF_ENERGY_TOL    ! Relative tolerance for energy stability
        real(dp) :: potential_tol = SCF_POTENTIAL_TOL ! Convergence tolerance for the potential residual
        real(dp) :: mixing_alpha = MIX_ALPHA       ! Mixing weight of the NEW potential. Used as-is when
                                                   ! use_adaptive_mixing = .false., and as the STARTING
                                                   ! value of the adaptive controller otherwise.
        logical :: verbose = .false.               ! Print convergence info
        logical :: store_history = .true.          ! Store convergence history
        logical :: use_adaptive_mixing = .true.    ! Use adaptive mixing (C++ behavior)
    end type scf_params_t

    !> @brief Results from SCF cycle
    type :: scf_results_t
        logical :: converged = .false.              ! Did SCF converge?
        integer :: n_iterations = 0                 ! Number of iterations performed
        real(dp) :: final_density_error = 0.0_dp    ! Final ||Δn||₂ (diagnostic only)
        real(dp) :: final_potential_residual = 0.0_dp ! Final ||V_calc - V_eff|| per site and spin
        real(dp) :: final_energy = 0.0_dp           ! Final total energy
        real(dp), allocatable :: density_up(:)      ! Converged spin-up density
        real(dp), allocatable :: density_down(:)    ! Converged spin-down density
        real(dp), allocatable :: eigvals(:)         ! Final eigenvalues (both spins)
        type(convergence_history_t) :: history      ! Convergence history
    end type scf_results_t

    !> Half width of the band around n = 1 in which a site counts as half filled
    real(dp), parameter, public :: HALF_FILLING_TOL = 1.0e-3_dp

    !> Number of consecutive sign alternations of ΔE required to call the energy
    !! "oscillating". One or two sign flips are common in healthy SCF runs; three
    !! in a row is a genuine period-2 oscillation.
    integer, parameter, public :: OSCILLATION_STREAK_MIN = 3

    public :: scf_params_t, scf_results_t
    public :: compute_total_energy
    public :: count_half_filled_sites, half_filling_warning_due
    public :: run_kohn_sham_scf_real
    public :: run_kohn_sham_scf_complex
    public :: init_scf_results, cleanup_scf_results

contains
    !> @brief Compute total Kohn-Sham energy
    !!
    !! E_tot = Σ_σ Σ_j ε_j,σ - U·Σ(n_up·n_down) - Σ(V_xc·n) + Σ(ε_xc)
    !!
    !! This formula matches the C++ original (lsdaks.cc lines 676-679).
    !!
    !! The eigenvalues ε_j already include V_eff = V_ext + U·n_other + V_xc.
    !! When summing eigenvalues, we get contributions from both spins:
    !!   - Spin up sees:   V_ext + U·n_down + V_xc_up
    !!   - Spin down sees: V_ext + U·n_up + V_xc_down
    !! This means the Hartree term U·n_other appears TWICE (once for each spin),
    !! giving 2·U·n_up·n_down instead of the correct U·n_up·n_down.
    !!
    !! The double-counting correction:
    !! - Band energy: Σε_j (includes 2·U·n_up·n_down)
    !! - Hartree correction: -U·Σ(n_up·n_down) (removes one copy)
    !! - V_xc correction: -Σ(V_xc·n) (removes V_xc from eigenvalues)
    !! - XC energy: +Σ(ε_xc) (adds true XC energy)
    !!
    !! Note: ε_xc from tables is the TOTAL XC energy at each site, not per particle!
    !!
    !! @param[in] eigvalues_up Spin-up eigenvalues (occupied ones)
    !! @param[in] eigvalues_down Spin-down eigenvalues (occupied ones)
    !! @param[in] n_up Number of spin-up electrons
    !! @param[in] n_down Number of spin-down electrons
    !! @param[in] density_up Spin-up density (length L)
    !! @param[in] density_down Spin-down density (length L)
    !! @param[in] V_ext External potential (length L)
    !! @param[in] xc_func XC functional object
    !! @param[in] U Hubbard interaction strength
    !! @param[in] L System size
    !! @param[out] total_energy Total energy
    !! @param[out] ierr Error code (0 = success)
    subroutine compute_total_energy(eigvals_up, eigvals_down, n_up, n_down, density_up, density_down, &
                                                            V_ext, xc_func, U, L, total_energy, ierr)
        real(dp), intent(in) :: eigvals_up(:), eigvals_down(:)
        integer, intent(in) :: n_up, n_down
        real(dp), intent(in) :: density_up(:), density_down(:), V_ext(:)
        type(xc_lsda_t), intent(in) :: xc_func
        real(dp), intent(in) :: U
        integer, intent(in) :: L
        real(dp), intent(out) :: total_energy
        integer, intent(out) :: ierr

        real(dp) :: E_band, E_hartree, E_xc_total, V_xc_correction
        real(dp) :: exc_val, V_xc_up, V_xc_down
        integer :: i

        ! 1. Band energy (sum of occupied eigenvalues)
        E_band = sum(eigvals_up(1:n_up)) + sum(eigvals_down(1:n_down))

        ! 2. Hartree energy: -U*Σ(n_up*n_down)
        !    Double-counting correction (matches C++ lsdaks.cc lines 676-679)
        !    Derivation: eigenvalues include V_eff = V_ext + U*n_other + V_xc
        !    When we sum v_eff*n for both spins, we get 2*U*n_up*n_down
        !    But we only want U*n_up*n_down, so we SUBTRACT U*n_up*n_down
        E_hartree = 0.0_dp
        do i = 1, L
            E_hartree = E_hartree + density_up(i) * density_down(i)
        end do
        E_hartree = -U * E_hartree  ! Negative sign!

        ! 3. Exchange-correlation energy and potential correction
        E_xc_total = 0.0_dp
        V_xc_correction = 0.0_dp

        do i = 1, L
            call get_exc(xc_func, density_up(i), density_down(i), exc_val, ierr)
            if (ierr /= ERROR_SUCCESS) then
                return
            end if

            call get_vxc(xc_func, density_up(i), density_down(i), &
                        V_xc_up, V_xc_down, ierr)

            if (ierr /= ERROR_SUCCESS) then
                return
            end if

            ! Note: exc_val is already the total XC energy at site i, not per particle!
            ! C++ just adds exc[i], not exc[i]*n_total[i]
            E_xc_total = E_xc_total + exc_val

            ! V_xc correction for double-counting
            V_xc_correction = V_xc_correction + V_xc_up * density_up(i) + V_xc_down * density_down(i)
        end do

        ! Total: E_band + E_hartree + E_xc - V_xc (matches C++ lsdaks.cc line 676-679)
        total_energy = E_band + E_hartree + E_xc_total - V_xc_correction

        ierr = ERROR_SUCCESS
    end subroutine compute_total_energy

    !> @brief Validate inputs for Kohn-Sham SCF cycle
    !!
    !! Performs comprehensive validation of system and SCF parameters before
    !! starting the self-consistent field iteration. Ensures physical consistency
    !! and prevents runtime errors from invalid inputs.
    !!
    !! Validation checks:
    !! - System size L > 0 (at least one lattice site)
    !! - Particle numbers per spin channel: 0 <= Nup <= L and 0 <= Ndown <= L
    !!   (Pauli exclusion: L orbitals per spin, hence N = Nup + Ndown <= 2L)
    !! - External potential size matches system size: size(V_ext) = L
    !! - SCF max_iter > 0 (at least one iteration allowed)
    !! - Mixing parameter 0 < mixing_alpha <= 1 (valid range for linear mixing)
    !! - potential_tol > 0 and energy_tol > 0. Both are convergence criteria
    !!   (residual_V < potential_tol AND |dE| < energy_tol*max(1,|E|)), and a
    !!   null or negative tolerance is unreachable: the cycle would burn the
    !!   whole iteration budget and report ERROR_CONVERGENCE_FAILED instead of
    !!   telling the caller that the request itself was impossible.
    !!
    !! @param[in]  params     System parameters (L, Nup, Ndown, bc, U, phase)
    !! @param[in]  scf_params SCF control parameters (max_iter, tolerances, mixing_alpha)
    !! @param[in]  V_ext      External potential array (length L)
    !! @param[out] ierr       Error code (ERROR_SUCCESS, ERROR_INVALID_INPUT, or ERROR_SIZE_MISMATCH)
    subroutine validate_kohn_sham_cycle_inputs(params, scf_params, V_ext, ierr)
        type(system_params_t), intent(in) :: params
        type(scf_params_t), intent(in) :: scf_params
        real(dp), intent(in) :: V_ext(:)
        integer, intent(out) :: ierr

        integer :: L, Nup, Ndown

        L = params%L
        Nup = params%Nup
        Ndown = params%Ndown

        ! Pauli exclusion per spin channel: each spin channel has exactly L
        ! orbitals, so 0 <= Nup <= L and 0 <= Ndown <= L. The total bound
        ! Nup + Ndown <= 2L follows from these two.
        if (L <= 0 .or. Nup < 0 .or. Ndown < 0 .or. Nup > L .or. Ndown > L) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (size(V_ext) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        if (scf_params%max_iter <= 0 .or. scf_params%mixing_alpha <= 0.0_dp .or. &
                                            scf_params%mixing_alpha > 1.0_dp) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        ! Both convergence tolerances must be strictly positive: residual_V and
        ! |E - E_prev| are non-negative, so a tolerance <= 0 can never be met.
        if (scf_params%potential_tol <= 0.0_dp) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (scf_params%energy_tol <= 0.0_dp) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        ierr = ERROR_SUCCESS
    end subroutine validate_kohn_sham_cycle_inputs

    !> @brief Count sites sitting at half filling
    !!
    !! A site counts as half filled when |n_up(i) + n_down(i) - 1| < HALF_FILLING_TOL.
    !! Those are exactly the sites that can cross the V_xc discontinuity at n = 1
    !! from one SCF iteration to the next (see `get_vxc`), flipping the local XC
    !! potential by 2|v_base(1, m)| and the total energy by the same amount.
    !!
    !! @param[in] density_up   Spin-up density (length >= L)
    !! @param[in] density_down Spin-down density (length >= L)
    !! @param[in] L            Number of lattice sites
    !! @return                 Number of sites within HALF_FILLING_TOL of n = 1
    pure function count_half_filled_sites(density_up, density_down, L) result(n_sites)
        real(dp), intent(in) :: density_up(:), density_down(:)
        integer, intent(in) :: L
        integer :: n_sites

        integer :: i

        n_sites = 0
        do i = 1, min(L, min(size(density_up), size(density_down)))
            if (abs(density_up(i) + density_down(i) - 1.0_dp) < HALF_FILLING_TOL) then
                n_sites = n_sites + 1
            end if
        end do
    end function count_half_filled_sites

    !> @brief Decide whether the half-filling / discontinuity warning applies
    !!
    !! The warning is meaningful only when BOTH conditions hold:
    !!   * at least one site sits within HALF_FILLING_TOL of n = 1, and
    !!   * the total energy is oscillating, i.e. the sign of ΔE has alternated
    !!     OSCILLATION_STREAK_MIN times in a row while |ΔE| is still above the
    !!     energy tolerance.
    !! Requiring both keeps the message off systems that merely sit near half
    !! filling and converge normally.
    !!
    !! @param[in] n_half_filled      Sites within HALF_FILLING_TOL of n = 1
    !! @param[in] alternation_streak Consecutive sign alternations of ΔE
    !! @param[in] delta_energy       Latest energy change E_iter - E_{iter-1}
    !! @param[in] total_energy       Latest total energy (relative scale)
    !! @param[in] energy_tol         Relative energy tolerance of the SCF cycle
    !! @return                       .true. if the warning should be printed
    pure function half_filling_warning_due(n_half_filled, alternation_streak, &
                                           delta_energy, total_energy, energy_tol) result(due)
        integer, intent(in) :: n_half_filled, alternation_streak
        real(dp), intent(in) :: delta_energy, total_energy, energy_tol
        logical :: due

        due = (n_half_filled > 0) .and. &
              (alternation_streak >= OSCILLATION_STREAK_MIN) .and. &
              (abs(delta_energy) > energy_tol * max(1.0_dp, abs(total_energy)))
    end function half_filling_warning_due

    !> @brief Run self-consistent Kohn-Sham cycle (real Hamiltonian)
    !!
    !! Iterates density → V_xc → H → diagonalize → new density until convergence.
    !! Mixing is applied to the POTENTIAL, never to the density.
    !!
    !! Convergence requires BOTH
    !!   * residual_V = ||V_calc - V_eff|| / sqrt(2L) < scf_params%potential_tol,
    !!     i.e. V_eff is a fixed point of the Kohn-Sham map, and
    !!   * |E - E_prev| < scf_params%energy_tol * max(1, |E|).
    !! ||Δn|| is reported but is not a criterion: it is proportional to the
    !! mixing weight and therefore vanishes whenever alpha does.
    !!
    !! @param[in] params System parameters (L, n_up, n_down, hopping, BC type)
    !! @param[in] scf_params SCF control parameters (max_iter, tolerances, mixing)
    !! @param[in] V_ext External potential V_ext(i) (length L)
    !! @param[in] xc_func XC functional object (already initialized with tables)
    !! @param[out] results SCF results (densities, eigenvalues, convergence info).
    !!                     This routine is the SOLE owner of the initialization of
    !!                     `results`: intent(out) resets every component on entry,
    !!                     the convergence history is allocated here (and only when
    !!                     scf_params%store_history is .true.), and on the rejection
    !!                     path of the validator `results` is returned untouched, with
    !!                     all its components deallocated. Callers must NOT pre-fill it
    !!                     with init_scf_results; that work would simply be discarded.
    !! @param[out] ierr Error code (0 = success)
    subroutine run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)
        type(system_params_t), intent(in) :: params
        type(scf_params_t), intent(in) :: scf_params
        real(dp), intent(in) :: V_ext(:)
        type(xc_lsda_t), intent(in) :: xc_func
        type(scf_results_t), intent(out) :: results
        integer, intent(out) :: ierr
        
        integer :: iter, i, L, Nup, Ndown
        real(dp) :: density_error, density_error_up, density_error_down, total_energy
        real(dp) :: residual_V, energy_prev
        logical :: is_converged, has_prev_energy, energy_is_stable
        real(dp) :: delta_energy, delta_energy_prev
        integer :: oscillation_streak, n_half_filled
        logical :: has_prev_delta, half_filling_warned
        type(adaptive_mix_t) :: mix_ctrl
        !> Mixing weight that actually produced the potential of the CURRENT
        !! iteration, and therefore its energy and its residuals. It is read
        !! from the controller at mixing time and never overwritten before the
        !! iteration is logged: the controller update at the end of the loop
        !! only decides the weight of the NEXT iteration.
        real(dp) :: alpha_used

        real(dp), allocatable :: n_up_in(:), n_down_in(:), n_up_out(:), n_down_out(:), V_xc_up(:), &
                            V_xc_down(:), H_up(:,:), H_down(:,:), eigvals_up(:), &
                            eigvals_down(:), eigvecs_up(:,:), eigvecs_down(:,:), delta_n_up(:), delta_n_down(:), &
                            V_eff_up(:), V_eff_down(:), V_eff_up_calc(:), V_eff_down_calc(:), &
                            V_zero(:)

        call validate_kohn_sham_cycle_inputs(params, scf_params, V_ext, ierr)

        if (ierr /= ERROR_SUCCESS) then
            return
        end if

        L = params%L
        Nup = params%Nup
        Ndown = params%Ndown

        allocate(n_up_in(L), n_down_in(L), n_up_out(L), n_down_out(L), V_xc_up(L), V_xc_down(L), &
                 H_up(L,L), H_down(L,L), eigvals_up(L), eigvals_down(L), &
                 eigvecs_up(L,L), eigvecs_down(L,L), delta_n_up(L), delta_n_down(L), &
                 V_eff_up(L), V_eff_down(L), V_eff_up_calc(L), V_eff_down_calc(L), V_zero(L))

        ! Initialize zero array for Hamiltonian builder
        V_zero(:) = 0.0_dp

        ! The history is allocated here and only here: this routine owns the
        ! initialization of `results` (see the intent(out) note above), and the
        ! arrays are only needed when the caller asked for them.
        if (scf_params%store_history) then
            call init_convergence_history(results%history, scf_params%max_iter, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for history init
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                           H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                           delta_n_down, V_eff_up, V_eff_down, V_eff_up_calc, V_eff_down_calc)
                return
            end if
        end if

        ! Initialize densities (uniform guess)
        n_up_in(:) = real(params%Nup, dp) / real(params%L, dp)
        n_down_in(:) = real(params%Ndown, dp) / real(params%L, dp)

        ! Initialize effective potentials from initial density guess
        ! This matches C++ initial_guess() function (lsdaks.cc lines 520-521)
        ! V_eff = V_ext + U*n_other + V_xc
        do i = 1, params%L
            ! Get V_xc from initial uniform density
            call get_vxc(xc_func, n_up_in(i), n_down_in(i), V_xc_up(i), V_xc_down(i), ierr)
            if (ierr /= ERROR_SUCCESS) then
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, V_eff_up, V_eff_down, V_eff_up_calc, V_eff_down_calc, V_zero)
                return
            end if

            ! Initialize V_eff = V_ext + U*n_other + V_xc (like C++)
            V_eff_up(i) = V_ext(i) + params%U * n_down_in(i) + V_xc_up(i)
            V_eff_down(i) = V_ext(i) + params%U * n_up_in(i) + V_xc_down(i)
        end do

        ! =================================
        ! Initialize adaptive mixing (if enabled)
        ! =================================
        if (scf_params%use_adaptive_mixing) then
            ! The adaptive controller starts from the user's mixing_alpha and
            ! retunes it from there. Seeding it with the hard-coded INITIAL_MIX
            ! instead would silently ignore scf_params%mixing_alpha whenever the
            ! controller is on, which is the default.
            call adaptive_mix_init(mix_ctrl, scf_params%energy_tol, &
                                   initial_alpha=scf_params%mixing_alpha)
            if (scf_params%verbose) then
                print '(A,F8.6)', "  Using adaptive mixing (C++ behavior), initial alpha = ", &
                    adaptive_mix_get_alpha(mix_ctrl)
            end if
        else
            if (scf_params%verbose) then
                print '(A,F6.4)', "  Using fixed mixing alpha = ", scf_params%mixing_alpha
            end if
        end if

        ! =================
        ! STEP 1: SCF loop
        ! =================

        residual_V = huge(1.0_dp)
        density_error = 0.0_dp
        density_error_up = 0.0_dp
        density_error_down = 0.0_dp
        total_energy = 0.0_dp
        energy_prev = 0.0_dp
        has_prev_energy = .false.
        delta_energy = 0.0_dp
        delta_energy_prev = 0.0_dp
        has_prev_delta = .false.
        oscillation_streak = 0
        n_half_filled = 0
        half_filling_warned = .false.

        do iter = 1, scf_params%max_iter
            ! -------------------------------------------------
            ! 1a. Compute V_xc from current densities
            ! -------------------------------------------------
            do i = 1, params%L
                call get_vxc(xc_func, n_up_in(i), n_down_in(i), V_xc_up(i), V_xc_down(i), ierr)
                if (ierr /= ERROR_SUCCESS) then
                    ! TODO: Proper error handling for V_xc calculation
                    deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, V_eff_up, V_eff_down, V_eff_up_calc, V_eff_down_calc, V_zero)
                    return
                end if
            end do

            ! -------------------------------------------------------------
            ! 1b. Calculate effective potentials V_eff = V_ext + U*n_other + V_xc
            ! (This is what C++ calls "v_ext[σ][j] + u*dens[other][j] + Vxc[σ][j]")
            ! -------------------------------------------------------------
            do i = 1, params%L
                ! V_eff_up_calc = V_ext + U*n_down + V_xc_up
                V_eff_up_calc(i) = V_ext(i) + params%U * n_down_in(i) + V_xc_up(i)

                ! V_eff_down_calc = V_ext + U*n_up + V_xc_down
                V_eff_down_calc(i) = V_ext(i) + params%U * n_up_in(i) + V_xc_down(i)
            end do

            ! -------------------------------------------------------------
            ! 1b'. True self-consistency residual of the effective potential
            !
            !   residual_V = sqrt( (||V_up_calc - V_up||^2 + ||V_down_calc - V_down||^2)
            !                      / (2*L) )
            !
            ! V_eff_up/V_eff_down still hold the potential that generated the
            ! current input density, so this is the distance between V_eff and
            ! its image under the Kohn-Sham map: it vanishes if and only if
            ! V_eff is a fixed point, and it does NOT scale with the mixing
            ! weight alpha (unlike ||Δn||, which is proportional to alpha because
            ! n_in is literally the previous n_out).
            !
            ! At iter = 1 the residual is identically zero by construction (V_eff
            ! is seeded from the same initial density), which is why convergence
            ! also requires a previous energy to compare against.
            ! -------------------------------------------------------------
            residual_V = sqrt((sum((V_eff_up_calc - V_eff_up)**2) + &
                               sum((V_eff_down_calc - V_eff_down)**2)) / real(2 * params%L, dp))

            ! -------------------------------------------------------------
            ! 1c. Mix potentials (C++ convention: Mix = weight of OLD)
            ! V_eff_new = Mix*V_eff_old + (1-Mix)*V_eff_calc
            ! -------------------------------------------------------------
            if (scf_params%use_adaptive_mixing) then
                ! Adaptive mixing uses the Mix parameter (updated each iteration)
                ! Convert to fortran alpha: alpha = 1 - Mix
                alpha_used = adaptive_mix_get_alpha(mix_ctrl)

                ! Apply mixing: V = (1-α)*V_old + α*V_calc = Mix*V_old + (1-Mix)*V_calc
                ! Since alpha = 1-Mix, we have: V = Mix*V_old + (1-Mix)*V_calc ✓
                do i = 1, params%L
                    V_eff_up(i) = (1.0_dp - alpha_used) * V_eff_up(i) + alpha_used * V_eff_up_calc(i)
                    V_eff_down(i) = (1.0_dp - alpha_used) * V_eff_down(i) + alpha_used * V_eff_down_calc(i)
                end do
            else
                ! Fixed mixing: alpha = weight of new
                alpha_used = scf_params%mixing_alpha
                do i = 1, params%L
                    V_eff_up(i) = (1.0_dp - alpha_used) * V_eff_up(i) + alpha_used * V_eff_up_calc(i)
                    V_eff_down(i) = (1.0_dp - alpha_used) * V_eff_down(i) + alpha_used * V_eff_down_calc(i)
                end do
            end if

            ! -------------------------------------
            ! 1d. Build Hamiltonians with mixed V_eff
            ! (C++ passes v_eff to hamiltonian_ks)
            ! -------------------------------------
            call build_hamiltonian(params%L, V_eff_up, V_zero, params%bc, params%phase, H_up, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for Hamiltonian Nup build
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, V_eff_up, V_eff_down, V_eff_up_calc, V_eff_down_calc, V_zero)
                return
            end if

            call build_hamiltonian(params%L, V_eff_down, V_zero, params%bc, params%phase, H_down, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for Hamiltonian Ndown build
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, V_eff_up, V_eff_down, V_eff_up_calc, V_eff_down_calc, V_zero)
                return
            end if

            ! ---------------------------------
            ! 1d. Diagonalize both Hamiltonians
            ! ---------------------------------
            call diagonalize_symmetric_real(H_up, params%L, eigvals_up, eigvecs_up, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for diagonalization Nup Hamiltonian
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down)
                return
            end if
            
            call diagonalize_symmetric_real(H_down, params%L, eigvals_down, eigvecs_down, ierr)
            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for diagonalization Ndown Hamiltonian
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down)
                return
            end if

            ! ----------------------------------------------
            ! 1e. Compute new densities from eigenvectors
            ! ----------------------------------------------
            call compute_density_spin(eigvecs_up, params%L, params%Nup, n_up_out, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for density calculation Nup
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down)
                return
            end if

            call compute_density_spin(eigvecs_down, params%L, params%Ndown, n_down_out, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for density calculation Ndown
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down)
                return
            end if
            
            ! ----------------------------------------------
            ! 1f. Compute density differences
            ! ----------------------------------------------
            call compute_density_difference(n_up_out, n_up_in, params%L, delta_n_up, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for Nup density difference
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down)
                return
            end if
            
            call compute_density_difference(n_down_out, n_down_in, params%L, delta_n_down, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for Ndown density difference
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down)
                return
            end if

            ! ---------------------------------------------------------------
            ! 1g. Density change, DIAGNOSTIC ONLY (see 1j for the criterion)
            !
            ! The two spin channels are measured separately and combined in
            ! quadrature. Summing Δn_up + Δn_down before taking the norm would
            ! cancel a pure polarisation oscillation (+d in the up channel, -d
            ! in the down channel) and report a spurious zero.
            ! ---------------------------------------------------------------
            call compute_density_norm(delta_n_up, params%L, L2, density_error_up, ierr)
            if (ierr == ERROR_SUCCESS) then
                call compute_density_norm(delta_n_down, params%L, L2, density_error_down, ierr)
            end if
            density_error = sqrt(density_error_up**2 + density_error_down**2)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for total density up norm
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down)
                return
            end if

                ! ------------------------
            ! 1h. Compute total energy
            ! ------------------------
            call compute_total_energy(eigvals_up, eigvals_down, params%Nup, params%Ndown, n_up_out, n_down_out, &
                                    V_ext, xc_func, params%U, params%L, total_energy, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for total energy calculation
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down)
                return
            end if

            ! ---------------------------------------------------------------
            ! 1h'. Half-filling / V_xc discontinuity diagnostic
            !
            ! V_xc is discontinuous at n = 1 (the BALDA Mott gap): a site that
            ! crosses that density between iterations flips its XC potential by
            ! 2|v_base(1, m)| and drags the total energy along, producing a
            ! period-2 oscillation that never settles. The condition below is
            ! deliberately conjunctive - half-filled sites AND an alternating,
            ! still-too-large ΔE - so that systems merely sitting near n = 1 and
            ! converging normally stay silent. The message is printed at most
            ! once per run.
            !
            ! The n = 1 discontinuity is NOT the only source of a persistent
            ! oscillation, and the warning says so: a partially filled,
            ! degenerate Fermi shell (e.g. L = 8, N_up = N_down = 4 with PBC)
            ! also fails to settle, because the density is built by filling the
            ! n_elec lowest eigenvectors with weight 1, without distributing the
            ! occupation over the degenerate level at the Fermi energy. That
            ! failure mode is independent of n = 1 (it also fires away from half
            ! filling) and smoothing V_xc does not cure it - it can make it
            ! worse. Fractional occupation is deferred to T7.
            ! ---------------------------------------------------------------
            if (has_prev_energy) then
                delta_energy = total_energy - energy_prev
                if (has_prev_delta) then
                    if (delta_energy * delta_energy_prev < 0.0_dp) then
                        oscillation_streak = oscillation_streak + 1
                    else
                        oscillation_streak = 0
                    end if
                end if
                delta_energy_prev = delta_energy
                has_prev_delta = .true.
            end if

            if (.not. half_filling_warned) then
                n_half_filled = count_half_filled_sites(n_up_out, n_down_out, params%L)
                if (half_filling_warning_due(n_half_filled, oscillation_streak, delta_energy, &
                                             total_energy, scf_params%energy_tol)) then
                    print '(A,I0,A)', "  WARNING: ", n_half_filled, &
                        " site(s) at half filling (|n - 1| < 1.0e-3); V_xc is discontinuous"
                    print '(A)', "           at n = 1 and the total energy may oscillate."
                    print '(A)', "           The oscillation may instead come from a partially filled," // &
                                 " degenerate Fermi shell,"
                    print '(A)', "           in which case smoothing V_xc does not help."
                    half_filling_warned = .true.
                end if
            end if

            ! ------------------
            ! 1i. Store history
            ! ------------------
            if (scf_params%store_history) then
                call update_convergence_history(iter, density_error, total_energy, results%history, ierr, &
                                                residual=residual_V)
            end if

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for history update
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down)
                return
            end if

            ! ------------------------------
            ! 1j. Update adaptive mixing (if enabled)
            !
            ! NOTE: this reads the weight for the NEXT iteration. The line
            ! logged below must report alpha_used - the weight that produced
            ! the residuals and the energy of THIS iteration - otherwise, on
            ! every up_mix/dw_mix transition, the log would attribute the
            ! numbers to a mixing weight that had nothing to do with them.
            ! ------------------------------
            if (scf_params%use_adaptive_mixing) then
                call adaptive_mix_update(mix_ctrl, total_energy)
            end if

            if (scf_params%verbose) then
                print '(A,I5,A,ES11.3,A,ES11.3,A,ES11.3,A,F16.8,A,F8.6)', "  Iter ", iter, &
                    "  |ΔV| = ", residual_V, "  |Δn↑| = ", density_error_up, &
                    "  |Δn↓| = ", density_error_down, "  E_tot = ", total_energy, &
                    "  α = ", alpha_used
            end if

            ! ---------------------------------------------------------------
            ! 1k. Convergence test
            !
            ! Self-consistency is declared only when the potential is a fixed
            ! point of the Kohn-Sham map (residual_V) AND the total energy has
            ! stopped moving in relative terms. ||Δn|| is NOT used: n_in is the
            ! previous n_out and the two come from potentials that differ by
            ! alpha*(V_calc - V_eff), so ||Δn|| is proportional to the mixing
            ! weight and shrinks whenever alpha shrinks, regardless of how far
            ! the cycle still is from self-consistency.
            !
            ! The energy comparison needs a previous iteration, which also rules
            ! out the trivially zero residual of iter = 1 (V_eff is seeded from
            ! the same density used to build V_eff_calc there).
            ! ---------------------------------------------------------------
            energy_is_stable = has_prev_energy .and. &
                abs(total_energy - energy_prev) < &
                    scf_params%energy_tol * max(1.0_dp, abs(total_energy))

            is_converged = (residual_V < scf_params%potential_tol) .and. energy_is_stable

            if (is_converged) then
                ! SUCCESS: Converged!
                results%converged = .true.
                results%n_iterations = iter
                results%final_density_error = density_error
                results%final_potential_residual = residual_V
                results%final_energy = total_energy

                ! Store final densities and eigenvalues
                allocate(results%density_up(params%L))
                allocate(results%density_down(params%L))
                allocate(results%eigvals(2*params%L))

                results%density_up = n_up_out
                results%density_down = n_down_out
                results%eigvals(1:params%L) = eigvals_up
                results%eigvals(params%L+1:2*params%L) = eigvals_down

                ierr = ERROR_SUCCESS
                return  ! return SCF LOOP
            end if

            energy_prev = total_energy
            has_prev_energy = .true.

            ! -----------------------------------------------------------------------
            ! 1l. Copy densities directly (NO MIXING!)
            ! -----------------------------------------------------------------------
            ! C++ code does: dens[0][i] = next_dens[0][i]  (line 695-696 in lsdaks.cc)
            ! Mixing is applied to POTENTIALS, not densities!
            n_up_in = n_up_out
            n_down_in = n_down_out
        end do

        ! ========================
        ! STEP 2: Did not converge
        ! ========================

        results%converged = .false.
        results%n_iterations = scf_params%max_iter
        results%final_density_error = density_error
        results%final_potential_residual = residual_V
        results%final_energy = total_energy

        ! The run is a failure, but the last densities and eigenvalues are still
        ! the best available estimate and the caller must be able to inspect them.
        allocate(results%density_up(params%L))
        allocate(results%density_down(params%L))
        allocate(results%eigvals(2*params%L))

        results%density_up = n_up_out
        results%density_down = n_down_out
        results%eigvals(1:params%L) = eigvals_up
        results%eigvals(params%L+1:2*params%L) = eigvals_down

        ierr = ERROR_CONVERGENCE_FAILED

        deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down)
    end subroutine run_kohn_sham_scf_real

    !> @brief Run self-consistent Kohn-Sham cycle (complex Hamiltonian)
    !!
    !! Iterates density → V_xc → H → diagonalize → new density until convergence.
    !! Mixing is applied to the POTENTIAL, never to the density.
    !!
    !! Convergence requires BOTH
    !!   * residual_V = ||V_calc - V_eff|| / sqrt(2L) < scf_params%potential_tol,
    !!     i.e. V_eff is a fixed point of the Kohn-Sham map, and
    !!   * |E - E_prev| < scf_params%energy_tol * max(1, |E|).
    !! ||Δn|| is reported but is not a criterion: it is proportional to the
    !! mixing weight and therefore vanishes whenever alpha does.
    !!
    !! @param[in] params System parameters (L, n_up, n_down, hopping, BC type)
    !! @param[in] scf_params SCF control parameters (max_iter, tolerances, mixing)
    !! @param[in] V_ext External potential V_ext(i) (length L)
    !! @param[in] xc_func XC functional object (already initialized with tables)
    !! @param[out] results SCF results (densities, eigenvalues, convergence info).
    !!                     This routine is the SOLE owner of the initialization of
    !!                     `results`: intent(out) resets every component on entry,
    !!                     the convergence history is allocated here (and only when
    !!                     scf_params%store_history is .true.), and on the rejection
    !!                     path of the validator `results` is returned untouched, with
    !!                     all its components deallocated. Callers must NOT pre-fill it
    !!                     with init_scf_results; that work would simply be discarded.
    !! @param[out] ierr Error code (0 = success)
    subroutine run_kohn_sham_scf_complex(params, scf_params, V_ext, xc_func, results, ierr)
        type(system_params_t), intent(in) :: params
        type(scf_params_t), intent(in) :: scf_params
        real(dp), intent(in) :: V_ext(:)
        type(xc_lsda_t), intent(in) :: xc_func
        type(scf_results_t), intent(out) :: results
        integer, intent(out) :: ierr

        integer :: iter, i, L, Nup, Ndown
        real(dp) :: density_error, density_error_up, density_error_down, total_energy
        real(dp) :: residual_V, energy_prev
        logical :: is_converged, has_prev_energy, energy_is_stable
        real(dp) :: delta_energy, delta_energy_prev
        integer :: oscillation_streak, n_half_filled
        logical :: has_prev_delta, half_filling_warned
        type(adaptive_mix_t) :: mix_ctrl
        !> Mixing weight that actually produced the potential of the CURRENT
        !! iteration, and therefore its energy and its residuals. It is read
        !! from the controller at mixing time and never overwritten before the
        !! iteration is logged: the controller update at the end of the loop
        !! only decides the weight of the NEXT iteration.
        real(dp) :: alpha_used

        real(dp), allocatable :: n_up_in(:), n_down_in(:), n_up_out(:), n_down_out(:), V_xc_up(:), &
                      V_xc_down(:), eigvals_up(:), eigvals_down(:), delta_n_up(:), delta_n_down(:), &
                      V_eff_up(:), V_eff_down(:), V_eff_up_calc(:), V_eff_down_calc(:), &
                      V_zero(:)

        complex(dp), allocatable :: H_up(:,:), H_down(:,:), eigvecs_up(:,:), eigvecs_down(:,:)

        call validate_kohn_sham_cycle_inputs(params, scf_params, V_ext, ierr)

        if (ierr /= ERROR_SUCCESS) then
            return
        end if

        L = params%L
        Nup = params%Nup
        Ndown = params%Ndown

        allocate(n_up_in(L), n_down_in(L), n_up_out(L), n_down_out(L), V_xc_up(L), V_xc_down(L), &
                 H_up(L,L), H_down(L,L), eigvals_up(L), eigvals_down(L), &
                 eigvecs_up(L,L), eigvecs_down(L,L), delta_n_up(L), delta_n_down(L), &
                 V_eff_up(L), V_eff_down(L), V_eff_up_calc(L), V_eff_down_calc(L), V_zero(L))

        ! Initialize zero array for Hamiltonian builder
        V_zero(:) = 0.0_dp

        ! The history is allocated here and only here: this routine owns the
        ! initialization of `results` (see the intent(out) note above), and the
        ! arrays are only needed when the caller asked for them.
        if (scf_params%store_history) then
            call init_convergence_history(results%history, scf_params%max_iter, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for history init
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                           H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                           delta_n_down, V_eff_up, V_eff_down, V_eff_up_calc, V_eff_down_calc, V_zero)
                return
            end if
        end if

        ! Initialize densities (uniform guess)
        n_up_in(:) = real(params%Nup, dp) / real(params%L, dp)
        n_down_in(:) = real(params%Ndown, dp) / real(params%L, dp)

        ! Initialize effective potentials from initial density guess
        ! This matches C++ initial_guess() function (lsdaks.cc lines 520-521)
        ! V_eff = V_ext + U*n_other + V_xc
        do i = 1, params%L
            ! Get V_xc from initial uniform density
            call get_vxc(xc_func, n_up_in(i), n_down_in(i), V_xc_up(i), V_xc_down(i), ierr)
            if (ierr /= ERROR_SUCCESS) then
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, V_eff_up, V_eff_down, V_eff_up_calc, V_eff_down_calc, V_zero)
                return
            end if

            ! Initialize V_eff = V_ext + U*n_other + V_xc (like C++)
            V_eff_up(i) = V_ext(i) + params%U * n_down_in(i) + V_xc_up(i)
            V_eff_down(i) = V_ext(i) + params%U * n_up_in(i) + V_xc_down(i)
        end do

        ! =================================
        ! Initialize adaptive mixing (if enabled)
        ! =================================
        if (scf_params%use_adaptive_mixing) then
            ! The adaptive controller starts from the user's mixing_alpha and
            ! retunes it from there. Seeding it with the hard-coded INITIAL_MIX
            ! instead would silently ignore scf_params%mixing_alpha whenever the
            ! controller is on, which is the default.
            call adaptive_mix_init(mix_ctrl, scf_params%energy_tol, &
                                   initial_alpha=scf_params%mixing_alpha)
            if (scf_params%verbose) then
                print '(A,F8.6)', "  Using adaptive mixing (C++ behavior), initial alpha = ", &
                    adaptive_mix_get_alpha(mix_ctrl)
            end if
        else
            if (scf_params%verbose) then
                print '(A,F6.4)', "  Using fixed mixing alpha = ", scf_params%mixing_alpha
            end if
        end if

        ! =================
        ! STEP 1: SCF loop
        ! =================

        residual_V = huge(1.0_dp)
        density_error = 0.0_dp
        density_error_up = 0.0_dp
        density_error_down = 0.0_dp
        total_energy = 0.0_dp
        energy_prev = 0.0_dp
        has_prev_energy = .false.
        delta_energy = 0.0_dp
        delta_energy_prev = 0.0_dp
        has_prev_delta = .false.
        oscillation_streak = 0
        n_half_filled = 0
        half_filling_warned = .false.

        do iter = 1, scf_params%max_iter
            ! -------------------------------------------------
            ! 1a. Compute V_xc from current densities
            ! -------------------------------------------------
            do i = 1, params%L
                call get_vxc(xc_func, n_up_in(i), n_down_in(i), V_xc_up(i), V_xc_down(i), ierr)
                if (ierr /= ERROR_SUCCESS) then
                    ! TODO: Proper error handling for V_xc calculation
                    deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, V_eff_up, V_eff_down, V_eff_up_calc, V_eff_down_calc, V_zero)
                    return
                end if
            end do

            ! -------------------------------------------------------------
            ! 1b. Calculate effective potentials V_eff = V_ext + U*n_other + V_xc
            ! (This is what C++ calls "v_ext[σ][j] + u*dens[other][j] + Vxc[σ][j]")
            ! -------------------------------------------------------------
            do i = 1, params%L
                ! V_eff_up_calc = V_ext + U*n_down + V_xc_up
                V_eff_up_calc(i) = V_ext(i) + params%U * n_down_in(i) + V_xc_up(i)

                ! V_eff_down_calc = V_ext + U*n_up + V_xc_down
                V_eff_down_calc(i) = V_ext(i) + params%U * n_up_in(i) + V_xc_down(i)
            end do

            ! -------------------------------------------------------------
            ! 1b'. True self-consistency residual of the effective potential
            !
            !   residual_V = sqrt( (||V_up_calc - V_up||^2 + ||V_down_calc - V_down||^2)
            !                      / (2*L) )
            !
            ! V_eff_up/V_eff_down still hold the potential that generated the
            ! current input density, so this is the distance between V_eff and
            ! its image under the Kohn-Sham map: it vanishes if and only if
            ! V_eff is a fixed point, and it does NOT scale with the mixing
            ! weight alpha (unlike ||Δn||, which is proportional to alpha because
            ! n_in is literally the previous n_out).
            !
            ! At iter = 1 the residual is identically zero by construction (V_eff
            ! is seeded from the same initial density), which is why convergence
            ! also requires a previous energy to compare against.
            ! -------------------------------------------------------------
            residual_V = sqrt((sum((V_eff_up_calc - V_eff_up)**2) + &
                               sum((V_eff_down_calc - V_eff_down)**2)) / real(2 * params%L, dp))

            ! -------------------------------------------------------------
            ! 1c. Mix potentials (C++ convention: Mix = weight of OLD)
            ! V_eff_new = Mix*V_eff_old + (1-Mix)*V_eff_calc
            ! -------------------------------------------------------------
            if (scf_params%use_adaptive_mixing) then
                ! Adaptive mixing uses the Mix parameter (updated each iteration)
                ! Convert to fortran alpha: alpha = 1 - Mix
                alpha_used = adaptive_mix_get_alpha(mix_ctrl)

                ! Apply mixing: V = (1-α)*V_old + α*V_calc = Mix*V_old + (1-Mix)*V_calc
                ! Since alpha = 1-Mix, we have: V = Mix*V_old + (1-Mix)*V_calc ✓
                do i = 1, params%L
                    V_eff_up(i) = (1.0_dp - alpha_used) * V_eff_up(i) + alpha_used * V_eff_up_calc(i)
                    V_eff_down(i) = (1.0_dp - alpha_used) * V_eff_down(i) + alpha_used * V_eff_down_calc(i)
                end do
            else
                ! Fixed mixing: alpha = weight of new
                alpha_used = scf_params%mixing_alpha
                do i = 1, params%L
                    V_eff_up(i) = (1.0_dp - alpha_used) * V_eff_up(i) + alpha_used * V_eff_up_calc(i)
                    V_eff_down(i) = (1.0_dp - alpha_used) * V_eff_down(i) + alpha_used * V_eff_down_calc(i)
                end do
            end if

            ! -------------------------------------
            ! 1d. Build Hamiltonians with mixed V_eff
            ! (C++ passes v_eff to hamiltonian_ks)
            ! -------------------------------------
            call build_hamiltonian_complex(params%L, V_eff_up, V_zero, params%bc, params%phase, H_up, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for Hamiltonian Nup build
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, V_eff_up, V_eff_down, V_eff_up_calc, V_eff_down_calc, V_zero)
                return
            end if

            call build_hamiltonian_complex(params%L, V_eff_down, V_zero, params%bc, params%phase, H_down, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for Hamiltonian Ndown build
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, V_eff_up, V_eff_down, V_eff_up_calc, V_eff_down_calc, V_zero)
                return
            end if

            ! ---------------------------------
            ! 1d. Diagonalize both Hamiltonians
            ! ---------------------------------
            call diagonalize_hermitian_complex(H_up, params%L, eigvals_up, eigvecs_up, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for diagonalization Nup Hamiltonian
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down)
                return
            end if
            
            call diagonalize_hermitian_complex(H_down, params%L, eigvals_down, eigvecs_down, ierr)
            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for diagonalization Ndown Hamiltonian
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down)
                return
            end if

            ! ----------------------------------------------
            ! 1e. Compute new densities from eigenvectors
            ! ----------------------------------------------
            call compute_density_spin(eigvecs_up, params%L, params%Nup, n_up_out, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for density calculation Nup
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down)
                return
            end if

            call compute_density_spin(eigvecs_down, params%L, params%Ndown, n_down_out, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for density calculation Ndown
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down)
                return
            end if
            
            ! ----------------------------------------------
            ! 1f. Compute density differences
            ! ----------------------------------------------
            call compute_density_difference(n_up_out, n_up_in, params%L, delta_n_up, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for Nup density difference
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down)
                return
            end if
            
            call compute_density_difference(n_down_out, n_down_in, params%L, delta_n_down, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for Ndown density difference
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down)
                return
            end if

            ! ---------------------------------------------------------------
            ! 1g. Density change, DIAGNOSTIC ONLY (see 1j for the criterion)
            !
            ! The two spin channels are measured separately and combined in
            ! quadrature. Summing Δn_up + Δn_down before taking the norm would
            ! cancel a pure polarisation oscillation (+d in the up channel, -d
            ! in the down channel) and report a spurious zero.
            ! ---------------------------------------------------------------
            call compute_density_norm(delta_n_up, params%L, L2, density_error_up, ierr)
            if (ierr == ERROR_SUCCESS) then
                call compute_density_norm(delta_n_down, params%L, L2, density_error_down, ierr)
            end if
            density_error = sqrt(density_error_up**2 + density_error_down**2)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for total density up norm
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down)
                return
            end if

                ! ------------------------
            ! 1h. Compute total energy
            ! ------------------------
            call compute_total_energy(eigvals_up, eigvals_down, params%Nup, params%Ndown, n_up_out, n_down_out, &
                                    V_ext, xc_func, params%U, params%L, total_energy, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for total energy calculation
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down)
                return
            end if

            ! ---------------------------------------------------------------
            ! 1h'. Half-filling / V_xc discontinuity diagnostic
            !
            ! V_xc is discontinuous at n = 1 (the BALDA Mott gap): a site that
            ! crosses that density between iterations flips its XC potential by
            ! 2|v_base(1, m)| and drags the total energy along, producing a
            ! period-2 oscillation that never settles. The condition below is
            ! deliberately conjunctive - half-filled sites AND an alternating,
            ! still-too-large ΔE - so that systems merely sitting near n = 1 and
            ! converging normally stay silent. The message is printed at most
            ! once per run.
            !
            ! The n = 1 discontinuity is NOT the only source of a persistent
            ! oscillation, and the warning says so: a partially filled,
            ! degenerate Fermi shell (e.g. L = 8, N_up = N_down = 4 with PBC)
            ! also fails to settle, because the density is built by filling the
            ! n_elec lowest eigenvectors with weight 1, without distributing the
            ! occupation over the degenerate level at the Fermi energy. That
            ! failure mode is independent of n = 1 (it also fires away from half
            ! filling) and smoothing V_xc does not cure it - it can make it
            ! worse. Fractional occupation is deferred to T7.
            ! ---------------------------------------------------------------
            if (has_prev_energy) then
                delta_energy = total_energy - energy_prev
                if (has_prev_delta) then
                    if (delta_energy * delta_energy_prev < 0.0_dp) then
                        oscillation_streak = oscillation_streak + 1
                    else
                        oscillation_streak = 0
                    end if
                end if
                delta_energy_prev = delta_energy
                has_prev_delta = .true.
            end if

            if (.not. half_filling_warned) then
                n_half_filled = count_half_filled_sites(n_up_out, n_down_out, params%L)
                if (half_filling_warning_due(n_half_filled, oscillation_streak, delta_energy, &
                                             total_energy, scf_params%energy_tol)) then
                    print '(A,I0,A)', "  WARNING: ", n_half_filled, &
                        " site(s) at half filling (|n - 1| < 1.0e-3); V_xc is discontinuous"
                    print '(A)', "           at n = 1 and the total energy may oscillate."
                    print '(A)', "           The oscillation may instead come from a partially filled," // &
                                 " degenerate Fermi shell,"
                    print '(A)', "           in which case smoothing V_xc does not help."
                    half_filling_warned = .true.
                end if
            end if

            ! ------------------
            ! 1i. Store history
            ! ------------------
            if (scf_params%store_history) then
                call update_convergence_history(iter, density_error, total_energy, results%history, ierr, &
                                                residual=residual_V)
            end if

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for history update
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down)
                return
            end if

            ! ------------------------------
            ! 1j. Update adaptive mixing (if enabled)
            !
            ! NOTE: this reads the weight for the NEXT iteration. The line
            ! logged below must report alpha_used - the weight that produced
            ! the residuals and the energy of THIS iteration - otherwise, on
            ! every up_mix/dw_mix transition, the log would attribute the
            ! numbers to a mixing weight that had nothing to do with them.
            ! ------------------------------
            if (scf_params%use_adaptive_mixing) then
                call adaptive_mix_update(mix_ctrl, total_energy)
            end if

            if (scf_params%verbose) then
                print '(A,I5,A,ES11.3,A,ES11.3,A,ES11.3,A,F16.8,A,F8.6)', "  Iter ", iter, &
                    "  |ΔV| = ", residual_V, "  |Δn↑| = ", density_error_up, &
                    "  |Δn↓| = ", density_error_down, "  E_tot = ", total_energy, &
                    "  α = ", alpha_used
            end if

            ! ---------------------------------------------------------------
            ! 1k. Convergence test
            !
            ! Self-consistency is declared only when the potential is a fixed
            ! point of the Kohn-Sham map (residual_V) AND the total energy has
            ! stopped moving in relative terms. ||Δn|| is NOT used: n_in is the
            ! previous n_out and the two come from potentials that differ by
            ! alpha*(V_calc - V_eff), so ||Δn|| is proportional to the mixing
            ! weight and shrinks whenever alpha shrinks, regardless of how far
            ! the cycle still is from self-consistency.
            !
            ! The energy comparison needs a previous iteration, which also rules
            ! out the trivially zero residual of iter = 1 (V_eff is seeded from
            ! the same density used to build V_eff_calc there).
            ! ---------------------------------------------------------------
            energy_is_stable = has_prev_energy .and. &
                abs(total_energy - energy_prev) < &
                    scf_params%energy_tol * max(1.0_dp, abs(total_energy))

            is_converged = (residual_V < scf_params%potential_tol) .and. energy_is_stable

            if (is_converged) then
                ! SUCCESS: Converged!
                results%converged = .true.
                results%n_iterations = iter
                results%final_density_error = density_error
                results%final_potential_residual = residual_V
                results%final_energy = total_energy

                ! Store final densities and eigenvalues
                allocate(results%density_up(params%L))
                allocate(results%density_down(params%L))
                allocate(results%eigvals(2*params%L))

                results%density_up = n_up_out
                results%density_down = n_down_out
                results%eigvals(1:params%L) = eigvals_up
                results%eigvals(params%L+1:2*params%L) = eigvals_down

                ierr = ERROR_SUCCESS
                return  ! return SCF LOOP
            end if

            energy_prev = total_energy
            has_prev_energy = .true.

            ! -----------------------------------------------------------------------
            ! 1l. Copy densities directly (NO MIXING!)
            ! -----------------------------------------------------------------------
            ! C++ code does: dens[0][i] = next_dens[0][i]  (line 695-696 in lsdaks.cc)
            ! Mixing is applied to POTENTIALS, not densities!
            n_up_in = n_up_out
            n_down_in = n_down_out
        end do

        ! ========================
        ! STEP 2: Did not converge
        ! ========================

        results%converged = .false.
        results%n_iterations = scf_params%max_iter
        results%final_density_error = density_error
        results%final_potential_residual = residual_V
        results%final_energy = total_energy

        ! The run is a failure, but the last densities and eigenvalues are still
        ! the best available estimate and the caller must be able to inspect them.
        allocate(results%density_up(params%L))
        allocate(results%density_down(params%L))
        allocate(results%eigvals(2*params%L))

        results%density_up = n_up_out
        results%density_down = n_down_out
        results%eigvals(1:params%L) = eigvals_up
        results%eigvals(params%L+1:2*params%L) = eigvals_down

        ierr = ERROR_CONVERGENCE_FAILED

        deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down)
    end subroutine run_kohn_sham_scf_complex

    !> @brief Initialize SCF results structure
    !!
    !! Helper for callers that assemble an scf_results_t outside the SCF cycle
    !! (the cycle itself owns the initialization of its own `results` and must
    !! NOT be pre-filled with this routine).
    !!
    !! The history allocation can fail - init_convergence_history rejects
    !! max_iter <= 0 with ERROR_INVALID_INPUT - and that failure is propagated.
    !! It used to be swallowed by an unconditional `ierr = ERROR_SUCCESS` at the
    !! end, so a caller asking for a history of zero iterations got a success
    !! code and an unallocated history.
    !!
    !! @param[out] results SCF results object
    !! @param[in] L System size
    !! @param[in] store_history Whether to allocate convergence history
    !! @param[in] max_iter Maximum iterations (for history size); must be > 0
    !!                     when store_history is .true.
    !! @param[out] ierr Error code (0 = success, ERROR_INVALID_INPUT if the
    !!                  history could not be allocated)
    subroutine init_scf_results(results, L, store_history, max_iter, ierr)
        type(scf_results_t), intent(out) :: results
        integer, intent(in) :: L, max_iter
        logical, intent(in) :: store_history
        integer, intent(out) :: ierr

        results%converged = .false.
        results%n_iterations = 0
        results%final_density_error = 0.0_dp
        results%final_potential_residual = 0.0_dp
        results%final_energy = 0.0_dp

        ierr = ERROR_SUCCESS

        if (store_history) then
            call init_convergence_history(results%history, max_iter, ierr)
            if (ierr /= ERROR_SUCCESS) return
        end if
    end subroutine init_scf_results

    !> @brief Deallocate SCF results
    !!
    !! @param[inout] results SCF results object
    !! @param[out] ierr Error code (0 = success)
    subroutine cleanup_scf_results(results, ierr)
        type(scf_results_t), intent(inout) :: results
        integer, intent(out) :: ierr

        if (allocated(results%density_up)) deallocate(results%density_up)
        if (allocated(results%density_down)) deallocate(results%density_down)
        if (allocated(results%eigvals)) deallocate(results%eigvals)

        call cleanup_convergence_history(results%history, ierr)
    end subroutine cleanup_scf_results

end module kohn_sham_cycle
