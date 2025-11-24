module kohn_sham_cycle
    use lsda_constants, only: dp, SCF_DENSITY_TOL, SCF_ENERGY_TOL, ITER_MAX, MIX_ALPHA
    use lsda_types, only: system_params_t
    use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT, &
                           ERROR_CONVERGENCE_FAILED, ERROR_SIZE_MISMATCH
    use boundary_conditions, only: apply_boundary_conditions, apply_boundary_conditions_complex
    use hamiltonian_builder, only: build_hamiltonian, build_hamiltonian_complex
    use lapack_wrapper, only: diagonalize_symmetric_real, diagonalize_hermitian_complex
    use density_calculator, only: compute_density_spin
    use xc_lsda, only: xc_lsda_t, get_vxc, get_exc
    use convergence_monitor, only: compute_density_difference, check_scf_convergence, &
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
        real(dp) :: density_tol = SCF_DENSITY_TOL  ! Convergence tolerance for density
        real(dp) :: energy_tol = SCF_ENERGY_TOL    ! Convergence tolerance for energy
        real(dp) :: mixing_alpha = MIX_ALPHA       ! Linear mixing parameter (used if use_adaptive_mixing = .false.)
        logical :: verbose = .false.               ! Print convergence info
        logical :: store_history = .true.          ! Store convergence history
        logical :: use_adaptive_mixing = .true.    ! Use adaptive mixing (C++ behavior)
    end type scf_params_t

    !> @brief Results from SCF cycle
    type :: scf_results_t
        logical :: converged                        ! Did SCF converge?
        integer :: n_iterations                     ! Number of iterations performed
        real(dp) :: final_density_error             ! Final ||Δn||₂
        real(dp) :: final_energy                    ! Final total energy
        real(dp), allocatable :: density_up(:)      ! Converged spin-up density
        real(dp), allocatable :: density_down(:)    ! Converged spin-down density
        real(dp), allocatable :: eigvals(:)         ! Final eigenvalues (both spins)
        type(convergence_history_t) :: history      ! Convergence history
    end type scf_results_t

    public :: scf_params_t, scf_results_t
    public :: compute_total_energy
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
    !! - Particle numbers Nup, Ndown >= 0 (non-negative)
    !! - Total particles N = Nup + Ndown <= L (Pauli exclusion: max 2 per site)
    !! - External potential size matches system size: size(V_ext) = L
    !! - SCF max_iter > 0 (at least one iteration allowed)
    !! - Mixing parameter 0 < mixing_alpha <= 1 (valid range for linear mixing)
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

        if (L <= 0 .or. Nup < 0 .or. Ndown < 0 .or. (Nup + Ndown) > 2*L) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (size(V_ext) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        if (scf_params%max_iter <= 0 .or. scf_params%mixing_alpha < 0.0_dp .or. &
                                            scf_params%mixing_alpha > 1.0_dp) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        ierr = ERROR_SUCCESS
    end subroutine validate_kohn_sham_cycle_inputs

    !> @brief Run self-consistent Kohn-Sham cycle (real Hamiltonian)
    !!
    !! Iterates density → V_xc → H → diagonalize → new density until convergence.
    !! Uses linear mixing to stabilize convergence.
    !!
    !! @param[in] params System parameters (L, n_up, n_down, hopping, BC type)
    !! @param[in] scf_params SCF control parameters (max_iter, tolerances, mixing)
    !! @param[in] V_ext External potential V_ext(i) (length L)
    !! @param[in] xc_func XC functional object (already initialized with tables)
    !! @param[out] results SCF results (densities, eigenvalues, convergence info)
    !! @param[out] ierr Error code (0 = success)
    subroutine run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)
        type(system_params_t), intent(in) :: params
        type(scf_params_t), intent(in) :: scf_params
        real(dp), intent(in) :: V_ext(:)
        type(xc_lsda_t), intent(in) :: xc_func
        type(scf_results_t), intent(out) :: results
        integer, intent(out) :: ierr
        
        integer :: iter, i, L, Nup, Ndown
        real(dp) :: density_error, total_energy
        logical :: is_converged
        type(adaptive_mix_t) :: mix_ctrl
        real(dp) :: current_alpha

        real(dp), allocatable :: n_up_in(:), n_down_in(:), n_up_out(:), n_down_out(:), V_xc_up(:), &
                            V_xc_down(:), H_up(:,:), H_down(:,:), eigvals_up(:), &
                            eigvals_down(:), eigvecs_up(:,:), eigvecs_down(:,:), delta_n_up(:), delta_n_down(:), &
                            delta_n_total(:), V_eff_up(:), V_eff_down(:), V_eff_up_calc(:), V_eff_down_calc(:), &
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
                 eigvecs_up(L,L), eigvecs_down(L,L), delta_n_up(L), delta_n_down(L), delta_n_total(L), &
                 V_eff_up(L), V_eff_down(L), V_eff_up_calc(L), V_eff_down_calc(L), V_zero(L))

        ! Initialize zero array for Hamiltonian builder
        V_zero(:) = 0.0_dp

        call init_convergence_history(results%history, scf_params%max_iter, ierr)

        if (ierr /= ERROR_SUCCESS) then
            ! TODO: Proper error handling for history init
            deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, delta_n_total, V_eff_up, V_eff_down, V_eff_up_calc, V_eff_down_calc)
            return
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
                       delta_n_down, delta_n_total, V_eff_up, V_eff_down, V_eff_up_calc, V_eff_down_calc, V_zero)
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
            call adaptive_mix_init(mix_ctrl, scf_params%energy_tol)
            if (scf_params%verbose) then
                print '(A)', "  Using adaptive mixing (C++ behavior)"
            end if
        else
            if (scf_params%verbose) then
                print '(A,F6.4)', "  Using fixed mixing alpha = ", scf_params%mixing_alpha
            end if
        end if

        ! =================
        ! STEP 1: SCF loop
        ! =================

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
                       delta_n_down, delta_n_total, V_eff_up, V_eff_down, V_eff_up_calc, V_eff_down_calc, V_zero)
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
            ! 1c. Mix potentials (C++ convention: Mix = weight of OLD)
            ! V_eff_new = Mix*V_eff_old + (1-Mix)*V_eff_calc
            ! -------------------------------------------------------------
            if (scf_params%use_adaptive_mixing) then
                ! Adaptive mixing uses the Mix parameter (updated each iteration)
                ! Convert to fortran alpha: alpha = 1 - Mix
                current_alpha = adaptive_mix_get_alpha(mix_ctrl)

                ! Apply mixing: V = (1-α)*V_old + α*V_calc = Mix*V_old + (1-Mix)*V_calc
                ! Since alpha = 1-Mix, we have: V = Mix*V_old + (1-Mix)*V_calc ✓
                do i = 1, params%L
                    V_eff_up(i) = (1.0_dp - current_alpha) * V_eff_up(i) + current_alpha * V_eff_up_calc(i)
                    V_eff_down(i) = (1.0_dp - current_alpha) * V_eff_down(i) + current_alpha * V_eff_down_calc(i)
                end do
            else
                ! Fixed mixing: alpha = weight of new
                current_alpha = scf_params%mixing_alpha
                do i = 1, params%L
                    V_eff_up(i) = (1.0_dp - current_alpha) * V_eff_up(i) + current_alpha * V_eff_up_calc(i)
                    V_eff_down(i) = (1.0_dp - current_alpha) * V_eff_down(i) + current_alpha * V_eff_down_calc(i)
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
                       delta_n_down, delta_n_total, V_eff_up, V_eff_down, V_eff_up_calc, V_eff_down_calc, V_zero)
                return
            end if

            call build_hamiltonian(params%L, V_eff_down, V_zero, params%bc, params%phase, H_down, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for Hamiltonian Ndown build
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, delta_n_total, V_eff_up, V_eff_down, V_eff_up_calc, V_eff_down_calc, V_zero)
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
                       delta_n_down, delta_n_total)
                return
            end if
            
            call diagonalize_symmetric_real(H_down, params%L, eigvals_down, eigvecs_down, ierr)
            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for diagonalization Ndown Hamiltonian
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, delta_n_total)
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
                       delta_n_down, delta_n_total)
                return
            end if

            call compute_density_spin(eigvecs_down, params%L, params%Ndown, n_down_out, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for density calculation Ndown
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, delta_n_total)
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
                       delta_n_down, delta_n_total)
                return
            end if
            
            call compute_density_difference(n_down_out, n_down_in, params%L, delta_n_down, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for Ndown density difference
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, delta_n_total)
                return
            end if

            ! ------------------------------
            ! 1g. Compute total density error
            ! ------------------------------
            delta_n_total = delta_n_up + delta_n_down
            call compute_density_norm(delta_n_total, params%L, L2, density_error, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for total density up norm
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, delta_n_total)
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
                       delta_n_down, delta_n_total)
                return
            end if

            ! ------------------
            ! 1i. Store history
            ! ------------------
            if (scf_params%store_history) then
                call update_convergence_history(iter, density_error, total_energy, results%history, ierr)
            end if

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for history update
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, delta_n_total)
                return
            end if

            ! ------------------------------
            ! 1j. Update adaptive mixing (if enabled)
            ! ------------------------------
            if (scf_params%use_adaptive_mixing) then
                call adaptive_mix_update(mix_ctrl, total_energy)
                current_alpha = adaptive_mix_get_alpha(mix_ctrl)

                if (scf_params%verbose) then
                    print '(A,I4,A,ES12.4,A,F16.8,A,F7.5)', "  Iter ", iter, "  |Δn| = ", density_error, &
                        "  E_tot = ", total_energy, "  α = ", current_alpha
                end if

                ! Check convergence: density-based (primary) OR energy-based (fallback)
                ! With potential mixing, density converges even if energy oscillates slightly
                call check_scf_convergence(delta_n_total, params%L, scf_params%density_tol, &
                                           is_converged, ierr)
                if (ierr /= ERROR_SUCCESS .or. .not. is_converged) then
                    ! Density not converged OR error - fall back to energy-based check
                    is_converged = mix_ctrl%converged
                end if
            else
                current_alpha = scf_params%mixing_alpha

                if (scf_params%verbose) then
                    print '(A,I4,A,ES12.4,A,F16.8)', "  Iter ", iter, "  |Δn| = ", density_error, &
                                                                        "  E_tot = ", total_energy
                end if

                ! 1k. Check convergence (fixed mixing)
                call check_scf_convergence(delta_n_total, params%L, scf_params%density_tol, &
                                           is_converged, ierr)
            end if

            ! Validate convergence check
            if (.not. scf_params%use_adaptive_mixing .and. ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for convergence check
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, delta_n_total)
                return
            end if

            if (is_converged) then
                ! SUCCESS: Converged!
                results%converged = .true.
                results%n_iterations = iter
                results%final_density_error = density_error
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
        results%final_energy = total_energy
        ierr = ERROR_CONVERGENCE_FAILED

        deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, delta_n_total)
    end subroutine run_kohn_sham_scf_real

    !> @brief Run self-consistent Kohn-Sham cycle (complex Hamiltonian)
    !!
    !! Iterates density → V_xc → H → diagonalize → new density until convergence.
    !! Uses linear mixing to stabilize convergence.
    !!
    !! @param[in] params System parameters (L, n_up, n_down, hopping, BC type)
    !! @param[in] scf_params SCF control parameters (max_iter, tolerances, mixing)
    !! @param[in] V_ext External potential V_ext(i) (length L)
    !! @param[in] xc_func XC functional object (already initialized with tables)
    !! @param[out] results SCF results (densities, eigenvalues, convergence info)
    !! @param[out] ierr Error code (0 = success)
    subroutine run_kohn_sham_scf_complex(params, scf_params, V_ext, xc_func, results, ierr)
        type(system_params_t), intent(in) :: params
        type(scf_params_t), intent(in) :: scf_params
        real(dp), intent(in) :: V_ext(:)
        type(xc_lsda_t), intent(in) :: xc_func
        type(scf_results_t), intent(out) :: results
        integer, intent(out) :: ierr

        integer :: iter, i, L, Nup, Ndown
        real(dp) :: density_error, total_energy
        logical :: is_converged
        type(adaptive_mix_t) :: mix_ctrl
        real(dp) :: current_alpha

        real(dp), allocatable :: n_up_in(:), n_down_in(:), n_up_out(:), n_down_out(:), V_xc_up(:), &
                      V_xc_down(:), eigvals_up(:), eigvals_down(:), delta_n_up(:), delta_n_down(:), &
                      delta_n_total(:), V_eff_up(:), V_eff_down(:), V_eff_up_calc(:), V_eff_down_calc(:), &
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
                 eigvecs_up(L,L), eigvecs_down(L,L), delta_n_up(L), delta_n_down(L), delta_n_total(L), &
                 V_eff_up(L), V_eff_down(L), V_eff_up_calc(L), V_eff_down_calc(L), V_zero(L))

        ! Initialize zero array for Hamiltonian builder
        V_zero(:) = 0.0_dp

        call init_convergence_history(results%history, scf_params%max_iter, ierr)

        if (ierr /= ERROR_SUCCESS) then
            ! TODO: Proper error handling for history init
            deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, delta_n_total, V_eff_up, V_eff_down, V_eff_up_calc, V_eff_down_calc, V_zero)
            return
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
                       delta_n_down, delta_n_total, V_eff_up, V_eff_down, V_eff_up_calc, V_eff_down_calc, V_zero)
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
            call adaptive_mix_init(mix_ctrl, scf_params%energy_tol)
            if (scf_params%verbose) then
                print '(A)', "  Using adaptive mixing (C++ behavior)"
            end if
        else
            if (scf_params%verbose) then
                print '(A,F6.4)', "  Using fixed mixing alpha = ", scf_params%mixing_alpha
            end if
        end if

        ! =================
        ! STEP 1: SCF loop
        ! =================

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
                       delta_n_down, delta_n_total, V_eff_up, V_eff_down, V_eff_up_calc, V_eff_down_calc, V_zero)
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
            ! 1c. Mix potentials (C++ convention: Mix = weight of OLD)
            ! V_eff_new = Mix*V_eff_old + (1-Mix)*V_eff_calc
            ! -------------------------------------------------------------
            if (scf_params%use_adaptive_mixing) then
                ! Adaptive mixing uses the Mix parameter (updated each iteration)
                ! Convert to fortran alpha: alpha = 1 - Mix
                current_alpha = adaptive_mix_get_alpha(mix_ctrl)

                ! Apply mixing: V = (1-α)*V_old + α*V_calc = Mix*V_old + (1-Mix)*V_calc
                ! Since alpha = 1-Mix, we have: V = Mix*V_old + (1-Mix)*V_calc ✓
                do i = 1, params%L
                    V_eff_up(i) = (1.0_dp - current_alpha) * V_eff_up(i) + current_alpha * V_eff_up_calc(i)
                    V_eff_down(i) = (1.0_dp - current_alpha) * V_eff_down(i) + current_alpha * V_eff_down_calc(i)
                end do
            else
                ! Fixed mixing: alpha = weight of new
                current_alpha = scf_params%mixing_alpha
                do i = 1, params%L
                    V_eff_up(i) = (1.0_dp - current_alpha) * V_eff_up(i) + current_alpha * V_eff_up_calc(i)
                    V_eff_down(i) = (1.0_dp - current_alpha) * V_eff_down(i) + current_alpha * V_eff_down_calc(i)
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
                       delta_n_down, delta_n_total, V_eff_up, V_eff_down, V_eff_up_calc, V_eff_down_calc, V_zero)
                return
            end if

            call build_hamiltonian_complex(params%L, V_eff_down, V_zero, params%bc, params%phase, H_down, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for Hamiltonian Ndown build
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, delta_n_total, V_eff_up, V_eff_down, V_eff_up_calc, V_eff_down_calc, V_zero)
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
                       delta_n_down, delta_n_total)
                return
            end if
            
            call diagonalize_hermitian_complex(H_down, params%L, eigvals_down, eigvecs_down, ierr)
            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for diagonalization Ndown Hamiltonian
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, delta_n_total)
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
                       delta_n_down, delta_n_total)
                return
            end if

            call compute_density_spin(eigvecs_down, params%L, params%Ndown, n_down_out, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for density calculation Ndown
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, delta_n_total)
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
                       delta_n_down, delta_n_total)
                return
            end if
            
            call compute_density_difference(n_down_out, n_down_in, params%L, delta_n_down, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for Ndown density difference
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, delta_n_total)
                return
            end if

            ! ------------------------------
            ! 1g. Compute total density error
            ! ------------------------------
            delta_n_total = delta_n_up + delta_n_down
            call compute_density_norm(delta_n_total, params%L, L2, density_error, ierr)

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for total density up norm
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, delta_n_total)
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
                       delta_n_down, delta_n_total)
                return
            end if

            ! ------------------
            ! 1i. Store history
            ! ------------------
            if (scf_params%store_history) then
                call update_convergence_history(iter, density_error, total_energy, results%history, ierr)
            end if

            if (ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for history update
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, delta_n_total)
                return
            end if

            ! ------------------------------
            ! 1j. Update adaptive mixing (if enabled)
            ! ------------------------------
            if (scf_params%use_adaptive_mixing) then
                call adaptive_mix_update(mix_ctrl, total_energy)
                current_alpha = adaptive_mix_get_alpha(mix_ctrl)

                if (scf_params%verbose) then
                    print '(A,I4,A,ES12.4,A,F16.8,A,F7.5)', "  Iter ", iter, "  |Δn| = ", density_error, &
                        "  E_tot = ", total_energy, "  α = ", current_alpha
                end if

                ! Check convergence: density-based (primary) OR energy-based (fallback)
                ! With potential mixing, density converges even if energy oscillates slightly
                call check_scf_convergence(delta_n_total, params%L, scf_params%density_tol, &
                                           is_converged, ierr)
                if (ierr /= ERROR_SUCCESS .or. .not. is_converged) then
                    ! Density not converged OR error - fall back to energy-based check
                    is_converged = mix_ctrl%converged
                end if
            else
                current_alpha = scf_params%mixing_alpha

                if (scf_params%verbose) then
                    print '(A,I4,A,ES12.4,A,F16.8)', "  Iter ", iter, "  |Δn| = ", density_error, &
                                                                        "  E_tot = ", total_energy
                end if

                ! 1k. Check convergence (fixed mixing)
                call check_scf_convergence(delta_n_total, params%L, scf_params%density_tol, &
                                           is_converged, ierr)
            end if

            ! Validate convergence check
            if (.not. scf_params%use_adaptive_mixing .and. ierr /= ERROR_SUCCESS) then
                ! TODO: Proper error handling for convergence check
                deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, delta_n_total)
                return
            end if

            if (is_converged) then
                ! SUCCESS: Converged!
                results%converged = .true.
                results%n_iterations = iter
                results%final_density_error = density_error
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
        results%final_energy = total_energy
        ierr = ERROR_CONVERGENCE_FAILED

        deallocate(n_up_in, n_down_in, n_up_out, n_down_out, V_xc_up, V_xc_down, &
                       H_up, H_down, eigvals_up, eigvals_down, eigvecs_up, eigvecs_down, delta_n_up, &
                       delta_n_down, delta_n_total)
    end subroutine run_kohn_sham_scf_complex

    !> @brief Initialize SCF results structure
    !!
    !! @param[out] results SCF results object
    !! @param[in] L System size
    !! @param[in] store_history Whether to allocate convergence history
    !! @param[in] max_iter Maximum iterations (for history size)
    !! @param[out] ierr Error code (0 = success)
    subroutine init_scf_results(results, L, store_history, max_iter, ierr)
        type(scf_results_t), intent(out) :: results
        integer, intent(in) :: L, max_iter
        logical, intent(in) :: store_history
        integer, intent(out) :: ierr

        results%converged = .false.
        results%n_iterations = 0
        results%final_density_error = 0.0_dp
        results%final_energy = 0.0_dp

        if (store_history) then
            call init_convergence_history(results%history, max_iter, ierr)
        end if

        ierr = ERROR_SUCCESS
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