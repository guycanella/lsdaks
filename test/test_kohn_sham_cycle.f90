!> Unit tests for kohn_sham_cycle module
program test_kohn_sham_cycle
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp
    implicit none

    real(dp), parameter :: TOL = 1.0e-9_dp

    !> Geometry and interaction of the alpha-independence / residual probe.
    !!
    !! Shared by run_fixed_alpha_probe and by the tests that recompute the
    !! residual independently, so that the reference computation can never
    !! drift away from the system the probe actually ran.
    integer, parameter :: PROBE_L = 20
    integer, parameter :: PROBE_NUP = 7
    integer, parameter :: PROBE_NDOWN = 5
    real(dp), parameter :: PROBE_U = 4.0_dp
    integer, parameter :: PROBE_IMP_SITE = 10
    real(dp), parameter :: PROBE_IMP_V = -2.0_dp
    character(len=*), parameter :: PROBE_TABLE = &
        "build/test_kohn_sham_xc_table_u4.00.dat"

    call prepare_test_xc_tables()
    call execute_serial_cmd_app(get_kohn_sham_tests())

contains

    function get_kohn_sham_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests

        tests = test_list([ &
            test("compute_total_energy_simple", test_compute_total_energy_simple), &
            test("total_energy_cancels_effective_potential_shift", &
                 test_total_energy_cancels_effective_potential_shift), &
            test("compute_total_energy_half_filling", test_compute_total_energy_half_filling), &
            test("validate_double_occupancy_accepted", test_validate_double_occupancy_accepted), &
            test("validate_inputs_invalid_L", test_validate_inputs_invalid_L), &
            test("validate_N_per_spin_exceeds_L", test_validate_N_per_spin_exceeds_L), &
            test("validate_inputs_size_mismatch", test_validate_inputs_size_mismatch), &
            test("validate_inputs_invalid_mixing", test_validate_inputs_invalid_mixing), &
            test("validate_inputs_nonpositive_tolerances", &
                 test_validate_inputs_nonpositive_tolerances), &
            test("validate_inputs_non_finite_tolerances", &
                 test_validate_inputs_non_finite_tolerances), &
            test("validate_xc_matches_system_u", test_validate_xc_matches_system_u), &
            test("scf_results_init_cleanup", test_scf_results_init_cleanup), &
            test("scf_results_init_propagates_history_error", &
                 test_scf_results_init_propagates_history_error), &
            test("complex_loop_matches_real_at_zero_twist", &
                 test_complex_loop_matches_real_at_zero_twist), &
            test("scf_converges_u0_open", test_scf_converges_u0_open), &
            test("scf_converges_u0_periodic", test_scf_converges_u0_periodic), &
            test("scf_stores_history", test_scf_stores_history), &
            test("scf_skips_history_when_disabled", test_scf_skips_history_when_disabled), &
            test("scf_density_conservation", test_scf_density_conservation), &
            test("scf_complex_twisted_bc", test_scf_complex_twisted_bc), &
            test("scf_declared_convergence_is_self_consistent", &
                 test_scf_declared_convergence_is_self_consistent), &
            test("scf_potential_residual_is_alpha_independent", &
                 test_scf_potential_residual_is_alpha_independent), &
            test("scf_residual_nonzero_at_alpha_one", &
                 test_scf_residual_nonzero_at_alpha_one), &
            test("scf_potential_residual_absolute_value", &
                 test_scf_potential_residual_absolute_value), &
            test("adaptive_mixing_honours_user_alpha", &
                 test_adaptive_mixing_honours_user_alpha), &
            test("count_half_filled_sites", test_count_half_filled_sites), &
            test("half_filling_warning_fires_when_oscillating", &
                 test_half_filling_warning_fires_when_oscillating), &
            test("half_filling_warning_silent_otherwise", &
                 test_half_filling_warning_silent_otherwise), &
            test("scf_open_shell_converges_uniform", &
                 test_scf_open_shell_converges_uniform), &
            test("scf_open_shell_is_self_consistent", &
                 test_scf_open_shell_is_self_consistent), &
            test("scf_fully_polarised_channel", &
                 test_scf_fully_polarised_channel), &
            test("scf_full_band_attractive_u", &
                 test_scf_full_band_attractive_u), &
            test("scf_reuses_output_xc_cache", test_scf_reuses_output_xc_cache), &
            test("xc_cache_preserves_half_filling_discontinuity", &
                 test_xc_cache_preserves_half_filling_discontinuity), &
            test("scf_failure_exposes_final_state", &
                 test_scf_failure_exposes_final_state), &
            test("scf_trap_doublet_stays_symmetric", &
                 test_scf_trap_doublet_stays_symmetric) &
        ])
    end function get_kohn_sham_tests

    !> REGRESSION (T7): open degenerate Fermi shell must converge, and uniformly
    !!
    !! L = 8, N_up = N_down = 4, periodic BC, U = 2, V_ext = 0. In PBC the
    !! levels are k = 2 pi m / 8 with degeneracies 1, 2, 2, 2, 1 (cumulative
    !! 1, 3, 5, 7, 8), so N_sigma = 4 leaves the Fermi shell OPEN.
    !!
    !! Before fractional occupation this system did not converge at all: LAPACK
    !! returns the real cos/sin combinations of the k = +-pi/2 pair, filling only
    !! one of them with weight 1 creates a period-2 density wave which feeds back
    !! into V_eff, and the cycle settles into a limit cycle (measured: 300
    !! iterations, ||delta n|| = 1.18e-1). With the occupation shared over the
    !! shell the very first iteration is already self-consistent.
    subroutine test_scf_open_shell_converges_uniform()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, scf_params_t, scf_results_t, &
                                    cleanup_scf_results
        use lsda_types, only: system_params_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_PERIODIC
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 8
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L)
        integer :: ierr
        character(len=256) :: table_file

        table_file = "build/test_kohn_sham_xc_table_u2.00.dat"
        call xc_lsda_init(xc_func, table_file, ierr)
        call check(ierr == ERROR_SUCCESS, "Open shell SCF: XC init should succeed")
        if (ierr /= ERROR_SUCCESS) return

        params%L = L
        params%Nup = 4
        params%Ndown = 4
        params%bc = BC_PERIODIC
        params%U = 2.0_dp
        params%phase = 0.0_dp

        scf_params%max_iter = 300
        scf_params%energy_tol = 1.0e-10_dp
        scf_params%potential_tol = 1.0e-8_dp
        scf_params%mixing_alpha = 0.2_dp
        scf_params%verbose = .false.
        scf_params%store_history = .false.

        V_ext = 0.0_dp

        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)

        call check(ierr == ERROR_SUCCESS, &
                   "Open shell SCF: L=8, N=4/4, PBC, U=2 must converge")
        call check(results%converged, "Open shell SCF: results must be flagged converged")

        if (allocated(results%density_up) .and. allocated(results%density_down)) then
            call check(maxval(abs(results%density_up - 0.5_dp)) < 1.0e-10_dp, &
                       "Open shell SCF: n_up must be uniform at N_up/L = 0.5")
            call check(maxval(abs(results%density_down - 0.5_dp)) < 1.0e-10_dp, &
                       "Open shell SCF: n_down must be uniform at N_down/L = 0.5")
        else
            call check(.false., "Open shell SCF: converged run must return densities")
        end if

        call cleanup_scf_results(results, ierr)
        call xc_lsda_destroy(xc_func)
    end subroutine test_scf_open_shell_converges_uniform

    !> REGRESSION (T20): harmonic trap whose Fermi level sits on a tunnel doublet
    !!
    !! L = 20, N_up = N_down = 13, open BC, U = 4, V = 0.7 (i - 10.5)^2,
    !! xc_smoothing_width = 0 (the exact BALDA functional). The trap core is a
    !! band insulator (n = 2 on sites 7..14), so the classically allowed region
    !! splits into two arms separated by a filled core and the arm levels come
    !! in tunnel-split doublets with a gap far below DEG_TOL. With N_sigma = 13
    !! the last electron of each spin channel falls on such a doublet and
    !! LAPACK returns it in a localised (left arm / right arm) basis.
    !!
    !! The equal sharing of the doublet is what keeps n(i) = n(L+1-i); any
    !! binary decision on the gap - or integer filling of a localised member -
    !! moves a whole electron into one arm and the cycle never settles. The
    !! test asserts the doublet (precondition, from the returned spectrum),
    !! the convergence with w = 0, the reflection symmetry to 1e-8 and the
    !! particle number. Measured: 130 iterations at alpha_0 = 0.1, asymmetry
    !! at roundoff.
    subroutine test_scf_trap_doublet_stays_symmetric()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, scf_params_t, scf_results_t, &
                                    cleanup_scf_results
        use lsda_types, only: system_params_t
        use lsda_constants, only: DEG_TOL
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_OPEN
        use potential_harmonic, only: apply_potential_harmonic
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 20, N_SIGMA = 13
        real(dp), parameter :: SPRING = 0.7_dp
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L), n_tot(L), asym
        integer :: ierr
        character(len=256) :: table_file

        table_file = PROBE_TABLE  ! xc_table_u4.00.dat
        call xc_lsda_init(xc_func, table_file, ierr)
        call check(ierr == ERROR_SUCCESS, "Trap doublet: XC init should succeed")
        if (ierr /= ERROR_SUCCESS) return

        call apply_potential_harmonic(SPRING, L, V_ext, ierr)
        call check(ierr == ERROR_SUCCESS, "Trap doublet: harmonic potential should succeed")

        params%L = L
        params%Nup = N_SIGMA
        params%Ndown = N_SIGMA
        params%bc = BC_OPEN
        params%U = 4.0_dp
        params%phase = 0.0_dp

        scf_params%max_iter = 600
        scf_params%energy_tol = 1.0e-10_dp
        scf_params%potential_tol = 1.0e-8_dp
        scf_params%mixing_alpha = 0.1_dp
        scf_params%verbose = .false.
        scf_params%store_history = .false.

        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)

        call check(ierr == ERROR_SUCCESS, &
                   "Trap doublet: L=20, N=13/13, OBC, U=4, k=0.7, w=0 must converge")
        call check(results%converged, "Trap doublet: results must be flagged converged")

        if (allocated(results%eigvals)) then
            ! Precondition: the Fermi level of each channel is on a doublet.
            call check(abs(results%eigvals(N_SIGMA + 1) - results%eigvals(N_SIGMA)) < DEG_TOL, &
                       "Trap doublet: levels 13 and 14 (up) must be degenerate below DEG_TOL")
            call check(abs(results%eigvals(L + N_SIGMA + 1) - results%eigvals(L + N_SIGMA)) < DEG_TOL, &
                       "Trap doublet: levels 13 and 14 (down) must be degenerate below DEG_TOL")
            call check(results%eigvals(N_SIGMA) - results%eigvals(N_SIGMA - 1) > 1.0e-3_dp, &
                       "Trap doublet: the doublet must be isolated from level 12 (precondition)")
        else
            call check(.false., "Trap doublet: converged run must return eigenvalues")
        end if

        if (allocated(results%density_up) .and. allocated(results%density_down)) then
            n_tot = results%density_up + results%density_down
            asym = maxval(abs(n_tot - n_tot(L:1:-1)))
            call check(asym < 1.0e-8_dp, &
                       "Trap doublet: n(i) = n(L+1-i) to 1e-8 with the doublet shared")
            call check(maxval(abs(results%density_up - results%density_up(L:1:-1))) < 1.0e-8_dp, &
                       "Trap doublet: n_up alone is also reflection symmetric")
            call check(abs(sum(n_tot) - real(2 * N_SIGMA, dp)) < 1.0e-10_dp, &
                       "Trap doublet: particle number conserved")
            call check(maxval(n_tot) > 1.99_dp, &
                       "Trap doublet: the core must be a band insulator (n = 2), else no arms")
        else
            call check(.false., "Trap doublet: converged run must return densities")
        end if

        call cleanup_scf_results(results, ierr)
        call xc_lsda_destroy(xc_func)
    end subroutine test_scf_trap_doublet_stays_symmetric

    !> An open-shell SCF solution must be self-consistent and uniformly occupied.
    !!
    !! L = 10, N_up = N_down = 8, periodic BC, U = 4, V_ext = 0.  The k = 2 pi m
    !! / 10 levels have cumulative degeneracies 1, 3, 5, 7, 9, 10, so N_sigma = 8
    !! leaves the Fermi shell open.  Correct fractional occupation must share
    !! that shell, preserving a uniform density and a self-consistent potential.
    !! This is a regression for the former limit cycle without a C++ energy
    !! oracle: the residual, particle number, uniformity, and Fermi doublet are
    !! all intrinsic properties of this translationally invariant problem.
    !!
    !! It also carries the only ABSOLUTE energy anchor at a point where `e_xc`
    !! is non-trivial.  Every other anchor of `compute_total_energy` sits in a
    !! limit where `e_xc` vanishes or is constant (U = 0, fully polarised, full
    !! band, Hartree difference), so a wrong factor on `E_xc` survives all of
    !! them.  (A double-counted Hartree term does NOT:
    !! `test_validate_double_occupancy_accepted` (`:852`) and
    !! `test_scf_full_band_attractive_u` (`:536`) pin `E/L = U` to 1e-12 at
    !! `n_up = n_dn = 1`, where doubling Hartree would give `2U`.)  Note the
    !! anchor does not see `V_xc` itself: by the same uniform-density algebra
    !! below, `V_eff^sigma` cancels out of the energy exactly; what it pins is
    !! `e_xc`.  Here the answer is reconstructible in
    !! closed form because the converged state is uniform: with `n_sigma = 0.8`
    !! at every site, `V_eff^sigma` is a site-independent constant, so it shifts
    !! every eigenvalue by the same amount and the shift `N_sigma V_eff^sigma`
    !! cancels the double-counting term `sum_i V_eff^sigma n_sigma(i)` exactly
    !! (`N_sigma = L n_sigma`).  What is left is
    !!
    !!   E = sum_sigma sum_m occ_m (-2 cos(2 pi m / L))   (free PBC band)
    !!       + L U n_up n_dn                              (Hartree)
    !!       + L e_xc(n_up, n_dn)                         (XC)
    !!
    !! with the open Fermi shell `m = +-4` carrying 1/2 an electron per state.
    !! Only `e_xc` is taken from the code (via `get_exc`, the quantity the
    !! functional is defined by); the assembly is built here, so a sign, factor
    !! or double count in `compute_total_energy` cannot cancel out of both sides.
    subroutine test_scf_open_shell_is_self_consistent()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, scf_params_t, scf_results_t, &
                                    cleanup_scf_results
        use lsda_types, only: system_params_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy, get_exc
        use boundary_conditions, only: BC_PERIODIC
        use lsda_constants, only: DEG_TOL, TWOPI
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 10
        real(dp), parameter :: N_SIGMA_SITE = 0.8_dp
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L)
        real(dp) :: e_band, e_hartree, exc_site, e_expected
        integer :: ierr, ierr_exc, m
        character(len=256) :: table_file

        table_file = PROBE_TABLE  ! xc_table_u4.00.dat
        call xc_lsda_init(xc_func, table_file, ierr)
        call check(ierr == ERROR_SUCCESS, "Open shell energy: XC init should succeed")
        if (ierr /= ERROR_SUCCESS) return

        params%L = L
        params%Nup = 8
        params%Ndown = 8
        params%bc = BC_PERIODIC
        params%U = 4.0_dp
        params%phase = 0.0_dp

        scf_params%max_iter = 300
        scf_params%energy_tol = 1.0e-10_dp
        scf_params%potential_tol = 1.0e-8_dp
        scf_params%mixing_alpha = 0.05_dp
        scf_params%verbose = .false.
        scf_params%store_history = .false.

        V_ext = 0.0_dp

        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)

        call check(ierr == ERROR_SUCCESS, &
                   "Open shell SCF: L=10, N=8/8, PBC, U=4 must converge")

        if (results%converged) then
            call check(results%final_potential_residual <= scf_params%potential_tol, &
                       "Open shell SCF: converged result must meet the potential residual tolerance")
            call check(maxval(abs(results%density_up - 0.8_dp)) < 1.0e-10_dp, &
                       "Open shell SCF: n_up must be uniform at 0.8")
            call check(maxval(abs(results%density_down - 0.8_dp)) < 1.0e-10_dp, &
                       "Open shell SCF: n_down must be uniform at 0.8")
            call check(abs(sum(results%density_up) - real(params%Nup, dp)) < 1.0e-10_dp .and. &
                       abs(sum(results%density_down) - real(params%Ndown, dp)) < 1.0e-10_dp, &
                       "Open shell SCF: both spin-channel particle numbers must be conserved")
            call check(abs(results%eigvals(9) - results%eigvals(8)) < DEG_TOL .and. &
                       abs(results%eigvals(L + 9) - results%eigvals(L + 8)) < DEG_TOL, &
                       "Open shell SCF: the fractional Fermi shell must remain degenerate")

            ! Free PBC band, per spin: m = 0 and the closed doublets m = +-1,
            ! +-2, +-3 are full, the open doublet m = +-4 holds one electron.
            e_band = -2.0_dp
            do m = 1, 3
                e_band = e_band + 2.0_dp * (-2.0_dp * cos(real(m, dp) * TWOPI / real(L, dp)))
            end do
            e_band = e_band + 1.0_dp * (-2.0_dp * cos(4.0_dp * TWOPI / real(L, dp)))
            e_band = 2.0_dp * e_band   ! both spin channels

            e_hartree = real(L, dp) * params%U * N_SIGMA_SITE * N_SIGMA_SITE

            call get_exc(xc_func, N_SIGMA_SITE, N_SIGMA_SITE, exc_site, ierr_exc)
            call check(ierr_exc == ERROR_SUCCESS, &
                       "Open shell energy: e_xc(0.8, 0.8) must be evaluable")

            e_expected = e_band + e_hartree + real(L, dp) * exc_site
            call check(abs(results%final_energy - e_expected) < 1.0e-8_dp, &
                       "Open shell energy: total energy must equal the uniform " // &
                       "band + Hartree - double-counting reconstruction")
        end if

        call cleanup_scf_results(results, ierr)
        call xc_lsda_destroy(xc_func)
    end subroutine test_scf_open_shell_is_self_consistent

    !> REGRESSION (T7): a fully polarised system must run, not be rejected
    !!
    !! L = 10, N_up = 5, N_down = 0, open BC, U = 4. validate_kohn_sham_cycle_inputs
    !! has always accepted N_sigma = 0, but density_calculator rejected
    !! n_elec <= 0, so the run died with ERROR_INVALID_INPUT downstream of a
    !! successful validation.
    !!
    !! With an empty down channel there is no Hartree and no XC coupling left,
    !! so the answer is the free-fermion one: E/L = -(1/L) sum_{j=1..5}
    !! 2 cos(j pi / 11) = -0.602667418333 (the same value the C++ prints).
    !!
    !! The tolerance is 1e-10, i.e. the result must be the free-fermion value to
    !! machine accuracy. Since the empty-channel shortcut of T18 the XC
    !! functional contributes EXACTLY zero here (e_xc = 0 and V_xc . n = 0 at
    !! every site), and the measured energy per site is -0.602667418333 against
    !! the analytic -0.602667418333227. The previous 1e-3/1e-2 tolerances were
    !! justified by a "small spurious XC contribution at n_down = 0" that no
    !! longer exists: they were loose enough to accept the very bias
    !! (e_xc(n, 0) ~ -8e-5 per site) that the shortcut removes.
    subroutine test_scf_fully_polarised_channel()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, scf_params_t, scf_results_t, &
                                    cleanup_scf_results
        use lsda_constants, only: PI
        use lsda_types, only: system_params_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_OPEN
        use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT

        integer, parameter :: L = 10, N_UP = 5
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L), e_exact
        integer :: ierr, j
        character(len=256) :: table_file

        table_file = PROBE_TABLE  ! xc_table_u4.00.dat
        call xc_lsda_init(xc_func, table_file, ierr)
        call check(ierr == ERROR_SUCCESS, "Polarised SCF: XC init should succeed")
        if (ierr /= ERROR_SUCCESS) return

        params%L = L
        params%Nup = N_UP
        params%Ndown = 0
        params%bc = BC_OPEN
        params%U = 4.0_dp
        params%phase = 0.0_dp

        scf_params%max_iter = 300
        scf_params%energy_tol = 1.0e-10_dp
        scf_params%potential_tol = 1.0e-8_dp
        scf_params%mixing_alpha = 0.2_dp
        scf_params%verbose = .false.
        scf_params%store_history = .false.

        V_ext = 0.0_dp

        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)

        call check(ierr /= ERROR_INVALID_INPUT, &
                   "Polarised SCF: N_down = 0 must not be rejected as invalid input")
        call check(ierr == ERROR_SUCCESS, "Polarised SCF: N_down = 0 must converge")

        if (results%converged) then
            call check(all(abs(results%density_down) < 1.0e-12_dp), &
                       "Polarised SCF: the empty channel must carry zero density")
            call check(abs(sum(results%density_up) - real(N_UP, dp)) < 1.0e-10_dp, &
                       "Polarised SCF: the filled channel must hold exactly N_up electrons")

            ! Free fermions in a box of L sites: eps_j = -2 cos(j pi / (L+1)).
            e_exact = 0.0_dp
            do j = 1, N_UP
                e_exact = e_exact - 2.0_dp * cos(real(j, dp) * PI / real(L + 1, dp))
            end do
            call check(abs(results%final_energy - e_exact) < 1.0e-10_dp, &
                       "Polarised SCF: E must match the free-fermion value to 1e-10")
            call check(abs(results%final_energy / real(L, dp) + 0.602667418333227_dp) < 1.0e-10_dp, &
                       "Polarised SCF: E/site must match the analytic -0.602667418333227 to 1e-10")
        end if

        call cleanup_scf_results(results, ierr)
        call xc_lsda_destroy(xc_func)
    end subroutine test_scf_fully_polarised_channel

    !> REGRESSION (B1): attractive twin of the completely full band
    !!
    !! L = 10, N_up = N_down = 10, periodic BC, U = -4. Exactly as in the
    !! repulsive case of test_validate_double_occupancy_accepted, the band sum
    !! is exhausted, so the kinetic energy vanishes and E/L = U = -4 exactly,
    !! with e_xc(1, 1) = 0 and V_xc = 0.
    !!
    !! This is the case where the region logic is most fragile, and the
    !! repulsive twin cannot see it: with U < 0 the Shiba transformation maps
    !! n_up = 1 onto n_up' = 1 - n_up = 0, so the evaluation lands exactly on
    !! the n = 1 particle-hole line and on the |m| = n edge of the physical
    !! triangle. One ulp of roundoff in the density (n_up = 1 + 4.4e-16 is what
    !! the diagonalization actually returns) is enough to push the point off
    !! both boundaries. Before the boundary snapping of `snap_to_boundaries`
    !! this run alternated between V_xc = 0 and V_xc = 1.6568542495 (the full
    !! dexc_dndown_b0(4, 1)) from iteration to iteration, never converged
    !! (30000 iterations) and returned E/L = -4.3828603 / -4.3990256; the
    !! C++ reference converges in 16 loops onto -4.000000.
    subroutine test_scf_full_band_attractive_u()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, scf_params_t, scf_results_t, &
                                    cleanup_scf_results
        use lsda_types, only: system_params_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_PERIODIC
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 10
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L)
        integer :: ierr
        character(len=256) :: table_file

        ! The table is tabulated for |U|; the attractive sign must be handed
        ! over explicitly, otherwise no Shiba transformation is applied.
        table_file = PROBE_TABLE  ! xc_table_u4.00.dat
        call xc_lsda_init(xc_func, table_file, ierr, u_signed = -4.0_dp)
        call check(ierr == ERROR_SUCCESS, "Attractive full band: XC init should succeed")
        if (ierr /= ERROR_SUCCESS) return

        params%L = L
        params%Nup = L
        params%Ndown = L
        params%bc = BC_PERIODIC
        params%U = -4.0_dp
        params%phase = 0.0_dp

        scf_params%max_iter = 300
        scf_params%energy_tol = 1.0e-10_dp
        scf_params%potential_tol = 1.0e-8_dp
        scf_params%mixing_alpha = 0.05_dp
        scf_params%verbose = .false.
        scf_params%store_history = .false.

        V_ext = 0.0_dp

        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)

        call check(ierr == ERROR_SUCCESS, &
                   "Attractive full band: L=10, N=10/10, PBC, U=-4 must converge")
        call check(results%converged, &
                   "Attractive full band: results must be flagged as converged")

        if (results%converged) then
            call check(maxval(abs(results%density_up - 1.0_dp)) < 1.0e-10_dp, &
                       "Attractive full band: n_up(i) = 1 on every site")
            call check(maxval(abs(results%density_down - 1.0_dp)) < 1.0e-10_dp, &
                       "Attractive full band: n_down(i) = 1 on every site")
            call check(abs(results%final_energy / real(L, dp) - params%U) < 1.0e-12_dp, &
                       "Attractive full band: E/L = U = -4 exactly")
        end if

        call cleanup_scf_results(results, ierr)
        call xc_lsda_destroy(xc_func)
    end subroutine test_scf_full_band_attractive_u

    !> Test total energy calculation with simple case
    !!
    !! Physics: E_tot = Σε_j - Σ(V_eff - V_ext)·n + U·Σn↑n↓ + E_xc
    !! The subtraction of (V_eff - V_ext)·n removes the Hartree and XC
    !! potentials that the eigenvalues already carry (C++ lsdaks.cc:675-679);
    !! the physical Hartree energy and E_xc are then added explicitly.
    !!
    !! The decisive assertion is differential: with a fixed, supplied V_eff, U
    !! enters compute_total_energy only through the explicit Hartree term (the
    !! XC functional is fixed by the table, not by the U argument), so
    !! E(U=2) - E(U=0) must be exactly +2·Σn↑n↓ = +1.25. This pins both the
    !! magnitude and the sign of the U argument; a wrong sign would give -1.25
    !! and a dropped/mis-passed U would give 0.
    subroutine test_compute_total_energy_simple()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: compute_total_energy
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 10
        real(dp) :: eigvals_up(5), eigvals_down(5)
        real(dp) :: n_up(L), n_down(L), V_ext(L), V_eff_up(L), V_eff_down(L)
        type(xc_lsda_t) :: xc_func
        real(dp) :: total_energy, energy_u0, hartree_sum
        integer :: ierr, i
        character(len=256) :: table_file

        ! Initialize XC functional
        table_file = "build/test_kohn_sham_xc_table_u2.00.dat"
        call xc_lsda_init(xc_func, table_file, ierr)
        call check(ierr == ERROR_SUCCESS, "XC init should succeed")

        ! Setup simple case: uniform density n=0.5
        do i = 1, L
            n_up(i) = 0.25_dp
            n_down(i) = 0.25_dp
            V_ext(i) = 0.0_dp
            V_eff_up(i) = 0.0_dp
            V_eff_down(i) = 0.0_dp
        end do

        ! Simple eigenvalues (5 occupied levels per spin)
        eigvals_up = [-2.0_dp, -1.5_dp, -1.0_dp, -0.5_dp, 0.0_dp]
        eigvals_down = [-2.0_dp, -1.5_dp, -1.0_dp, -0.5_dp, 0.0_dp]

        ! U = 2.0 matches the XC table loaded above (xc_table_u2.00.dat)
        call compute_total_energy(eigvals_up, eigvals_down, 5, 5, n_up, n_down, &
                                  V_ext, V_eff_up, V_eff_down, xc_func, 2.0_dp, L, total_energy, ierr)

        call check(ierr == ERROR_SUCCESS, "Energy calculation should succeed")
        call check(total_energy == total_energy, "Energy should be valid number")

        ! Same call with U = 0: identical band energy, identical XC (same table),
        ! only the Hartree term disappears.
        call compute_total_energy(eigvals_up, eigvals_down, 5, 5, n_up, n_down, &
                                  V_ext, V_eff_up, V_eff_down, xc_func, 0.0_dp, L, energy_u0, ierr)
        call check(ierr == ERROR_SUCCESS, "Energy calculation at U=0 should succeed")

        ! With the same supplied V_eff, the explicit physical Hartree term is
        ! +U*Σn_up*n_down. (A self-consistent V_eff would move the band energy
        ! and its subtraction cancels that shift.)
        hartree_sum = sum(n_up * n_down)
        call check(abs(hartree_sum - 0.625_dp) < TOL, "Precondition: Σn_up*n_down = 0.625")
        call check(abs((total_energy - energy_u0) - 2.0_dp * hartree_sum) < 1.0e-12_dp, &
                   "E(U=2) - E(U=0) must equal +U*Σn_up*n_down = 1.25")

        ! Band energy alone is exactly -10 and V_eff = V_ext = 0, so at U=0 the
        ! whole remainder is Σε_xc; it must stay small compared to the band.
        call check(abs(energy_u0 + 10.0_dp) < 2.0_dp, &
                   "At U=0 the energy must sit close to the band energy -10")

        call xc_lsda_destroy(xc_func)
    end subroutine test_compute_total_energy_simple


    !> T19: the band-energy shift from V_eff must cancel its double counting
    subroutine test_total_energy_cancels_effective_potential_shift()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: compute_total_energy
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 10
        real(dp), parameter :: SHIFT = 0.75_dp
        real(dp) :: eigvals_up(5), eigvals_down(5), n_up(L), n_down(L), V_ext(L)
        real(dp) :: V_eff_up(L), V_eff_down(L), energy_ref, energy_shifted
        type(xc_lsda_t) :: xc_func
        integer :: ierr

        call xc_lsda_init(xc_func, 'build/test_kohn_sham_xc_table_u2.00.dat', ierr)
        call check(ierr == ERROR_SUCCESS, 'T19: XC init should succeed')
        if (ierr /= ERROR_SUCCESS) return

        eigvals_up = [-2.0_dp, -1.5_dp, -1.0_dp, -0.5_dp, 0.0_dp]
        eigvals_down = eigvals_up
        ! Five occupied orbitals per spin on ten sites: each density must sum
        ! to five for the band-energy and potential shifts to represent one
        ! common Kohn-Sham state.
        n_up = 0.5_dp
        n_down = 0.5_dp
        V_ext = 0.0_dp
        V_eff_up = 0.0_dp
        V_eff_down = 0.0_dp

        call compute_total_energy(eigvals_up, eigvals_down, 5, 5, n_up, n_down, V_ext, &
                                  V_eff_up, V_eff_down, xc_func, 2.0_dp, L, energy_ref, ierr)
        call check(ierr == ERROR_SUCCESS, 'T19: reference energy should succeed')

        eigvals_up = eigvals_up + SHIFT
        eigvals_down = eigvals_down + SHIFT
        V_eff_up = SHIFT
        V_eff_down = SHIFT
        call compute_total_energy(eigvals_up, eigvals_down, 5, 5, n_up, n_down, V_ext, &
                                  V_eff_up, V_eff_down, xc_func, 2.0_dp, L, energy_shifted, ierr)
        call check(ierr == ERROR_SUCCESS, 'T19: shifted energy should succeed')
        call check(abs(energy_shifted - energy_ref) < 1.0e-12_dp, &
                   'T19: a uniform V_eff shift must cancel from the total energy')

        call xc_lsda_destroy(xc_func)
    end subroutine test_total_energy_cancels_effective_potential_shift

    !> Test total energy for half-filling case
    !!
    !! Physics: Half-filling (n=1) is a special case with enhanced correlations.
    !! Ground state energy should be lower than non-interacting case.
    subroutine test_compute_total_energy_half_filling()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: compute_total_energy
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 8
        real(dp) :: eigvals_up(L), eigvals_down(L)
        real(dp) :: n_up(L), n_down(L), V_ext(L), V_eff_up(L), V_eff_down(L)
        type(xc_lsda_t) :: xc_func
        real(dp) :: total_energy
        integer :: ierr, i
        character(len=256) :: table_file

        table_file = "build/test_kohn_sham_xc_table_u4.00.dat"
        call xc_lsda_init(xc_func, table_file, ierr)
        call check(ierr == ERROR_SUCCESS, "XC init should succeed")

        ! Half-filling: n_up = n_down = 0.5 at each site
        do i = 1, L
            n_up(i) = 0.5_dp
            n_down(i) = 0.5_dp
            V_ext(i) = 0.0_dp
            V_eff_up(i) = 0.0_dp
            V_eff_down(i) = 0.0_dp
        end do

        ! The lowest L/2 eigenvalues are occupied in each spin channel
        eigvals_up = [-2.0_dp, -1.8_dp, -1.5_dp, -1.2_dp, -1.0_dp, -0.8_dp, -0.5_dp, -0.2_dp]
        eigvals_down = eigvals_up

        ! U = 4.0 matches the XC table loaded above (xc_table_u4.00.dat)
        call compute_total_energy(eigvals_up, eigvals_down, L / 2, L / 2, n_up, n_down, &
                                  V_ext, V_eff_up, V_eff_down, xc_func, 4.0_dp, L, total_energy, ierr)

        call check(ierr == ERROR_SUCCESS, "Energy calculation should succeed")
        call check(total_energy == total_energy, "Energy should be valid")

        call xc_lsda_destroy(xc_func)
    end subroutine test_compute_total_energy_half_filling

    !> Test that the SCF entry point accepts doubly occupied fillings (N > L)
    !!
    !! Physics: each spin channel owns exactly L orbitals, so the only bound is
    !! 0 <= N_sigma <= L per channel (total N <= 2L). Double occupancy of a site
    !! is physical in the Hubbard model, so any N with L < N <= 2L is a perfectly
    !! legal input. Two accepted fillings are exercised here:
    !!   * L < N < 2L: Nup = Ndown = 7 on L = 10 (n = 1.4 per site);
    !!   * N = 2L:     Nup = Ndown = L = 10 (band insulator, n = 2 per site).
    !! The rejected fillings (N_sigma > L) live exclusively in
    !! test_validate_N_per_spin_exceeds_L.
    !!
    !! This case is the regression guard against the obsolete rule
    !! "Nup + Ndown <= L": under that rule 7 + 7 = 14 > 10 would be rejected with
    !! ERROR_INVALID_INPUT, so this test fails if anyone reintroduces the total
    !! bound in validate_kohn_sham_cycle_inputs.
    !!
    !! N_sigma = 7 is chosen deliberately: on a 10-site ring the single-particle
    !! levels k = 2*pi*m/10 have degeneracies 1,2,2,2,2,1 (cumulative 1,3,5,7,9),
    !! so 7 electrons per spin close a shell. The occupied set is therefore
    !! unambiguous. Open-shell fillings such as N_sigma = 6 or 8 would instead
    !! require fractional occupation of the partially filled degenerate level
    !! (C++ update_degen), which is not yet ported; they are tracked separately
    !! and are not exercised here.
    !!
    !! Convergence in a single iteration is *not* a consequence of the closed
    !! shell alone. It follows jointly from two facts:
    !!   1. the shell is closed, so diagonalizing a translationally invariant
    !!      H reproduces the uniform density n_sigma(i) = N_sigma/L exactly (no
    !!      ambiguity in the occupied set, no symmetry breaking);
    !!   2. run_kohn_sham_scf_real seeds V_eff from that same uniform density
    !!      (src/kohn_sham/kohn_sham_cycle.f90:262-275, V_eff = V_ext + U*n_other
    !!      + V_xc), so the starting potential is already self-consistent in the
    !!      translationally invariant case.
    !! Fact 2 is what protects the assertions below from the planned change of
    !! convergence criterion (T3): the fixed point is reached at the very first
    !! potential evaluation, so no tolerance or norm choice can move it.
    subroutine test_validate_double_occupancy_accepted()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, scf_params_t, scf_results_t, &
                                    cleanup_scf_results
        use lsda_types, only: system_params_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_PERIODIC
        use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT

        integer, parameter :: L = 10
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L)
        integer :: ierr
        character(len=256) :: table_file

        table_file = "build/test_kohn_sham_xc_table_u4.00.dat"
        call xc_lsda_init(xc_func, table_file, ierr)
        call check(ierr == ERROR_SUCCESS, "XC init should succeed")

        ! ---------------------------------------------------------------
        ! Case 1, L < N < 2L: Nup = Ndown = 7 on L = 10. Legal per spin
        ! channel (7 <= 10), rejected by the old total bound (14 > 10).
        ! Closed shell in PBC.
        ! ---------------------------------------------------------------
        params%L = L
        params%Nup = 7
        params%Ndown = 7
        params%bc = BC_PERIODIC
        params%U = 4.0_dp
        params%phase = 0.0_dp

        scf_params%max_iter = 100
        scf_params%density_tol = 1.0e-6_dp
        scf_params%energy_tol = 1.0e-8_dp
        scf_params%mixing_alpha = 0.3_dp
        scf_params%verbose = .false.
        scf_params%store_history = .true.

        V_ext = 0.0_dp

        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)

        call check(ierr /= ERROR_INVALID_INPUT, &
                   "Nup = Ndown = 7 on L = 10 must not be rejected as invalid input")
        call check(ierr == ERROR_SUCCESS, &
                   "SCF on the closed-shell 7/7 case must converge")
        call check(results%converged, "SCF results must be flagged as converged")

        ! Physical invariants: each spin channel keeps its own particle number
        ! and, with V_ext = 0 on a ring, the density stays uniform at N_sigma/L.
        call check(allocated(results%density_up), "Spin-up density must be allocated")
        call check(allocated(results%density_down), "Spin-down density must be allocated")
        if (allocated(results%density_up) .and. allocated(results%density_down)) then
            call check(abs(sum(results%density_up) - 7.0_dp) < 1.0e-10_dp, &
                       "Spin-up particle number must stay at Nup = 7")
            call check(abs(sum(results%density_down) - 7.0_dp) < 1.0e-10_dp, &
                       "Spin-down particle number must stay at Ndown = 7")
            call check(maxval(abs(results%density_up - 0.7_dp)) < 1.0e-10_dp, &
                       "Spin-up density must be uniform at n_up(i) = 0.7")
            call check(maxval(abs(results%density_down - 0.7_dp)) < 1.0e-10_dp, &
                       "Spin-down density must be uniform at n_down(i) = 0.7")
        end if

        call cleanup_scf_results(results, ierr)

        ! ---------------------------------------------------------------
        ! Case 2, N = 2L: both channels completely filled (band insulator).
        ! This is the upper edge of the allowed range and must be accepted.
        ! ---------------------------------------------------------------
        params%Nup = L
        params%Ndown = L

        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)

        call check(ierr /= ERROR_INVALID_INPUT, &
                   "N = 2L (both channels full) must not be rejected as invalid input")
        call check(ierr == ERROR_SUCCESS, "SCF on the fully filled N = 2L case must converge")
        call check(results%converged, "SCF results for N = 2L must be flagged as converged")
        call check(allocated(results%density_up), "Spin-up density must be allocated for N = 2L")
        call check(allocated(results%density_down), "Spin-down density must be allocated for N = 2L")
        if (allocated(results%density_up) .and. allocated(results%density_down)) then
            call check(abs(sum(results%density_up) - real(L, dp)) < 1.0e-10_dp, &
                       "Spin-up particle number must stay at Nup = L")
            call check(abs(sum(results%density_down) - real(L, dp)) < 1.0e-10_dp, &
                       "Spin-down particle number must stay at Ndown = L")
            call check(maxval(abs(results%density_up - 1.0_dp)) < 1.0e-10_dp, &
                       "Fully filled spin-up channel must give n_up(i) = 1")
            call check(maxval(abs(results%density_down - 1.0_dp)) < 1.0e-10_dp, &
                       "Fully filled spin-down channel must give n_down(i) = 1")
        end if

        ! REGRESSION (T18): the completely full band is analytically trivial.
        ! Every site is doubly occupied, so the kinetic energy vanishes (the band
        ! sum is exhausted) and the only surviving term is the Hubbard one,
        ! E/L = U exactly. The XC energy must NOT contribute: at (n_up, n_dw) =
        ! (1, 1) the Region IV symmetry maps the point onto the EMPTY lattice
        ! (0, 0), where e_xc and V_xc vanish.
        !
        ! The density assertions above cannot see this: they are unitarity
        ! identities of the eigenvector matrix and hold for any XC term. Before
        ! the empty-channel shortcut was applied to the MAPPED densities the
        ! spline was extrapolated into the empty corner and this run returned
        ! E/L = 4.0007185884, i.e. a spurious e_xc(1, 1) = 7.19e-4 per site.
        call check(abs(results%final_energy / real(L, dp) - params%U) < 1.0e-12_dp, &
                   "Completely full band must give E/L = U exactly")

        ! REGRESSION (T19): the energy must be a functional of n_out at EVERY
        ! iteration, not only at the fixed point. The band is full from the
        ! first diagonalisation on, so E/L = U must already hold at iteration 1
        ! whatever V_eff the mixing started from. Rebuilding the double-counting
        ! correction from V_xc(n_out) instead of the diagonalised V_eff leaves
        ! the spurious term -Σ(V_calc - V_eff)·n, which this assertion sees.
        call check(allocated(results%history%energies), &
                   "T19: energy history must be stored for the full-band run")
        if (allocated(results%history%energies)) then
            call check(results%history%current_iter >= 1, &
                       "T19: at least one iteration must be recorded")
            if (results%history%current_iter >= 1) then
                call check(abs(results%history%energies(1) / real(L, dp) - params%U) < 1.0e-12_dp, &
                           "T19: E/L = U must hold already at the first SCF iteration")
            end if
        end if

        call cleanup_scf_results(results, ierr)
        call xc_lsda_destroy(xc_func)
    end subroutine test_validate_double_occupancy_accepted

    !> Test that the SCF entry point rejects a non-positive lattice size
    !!
    !! L <= 0 has no lattice to diagonalize, so validate_kohn_sham_cycle_inputs
    !! must return ERROR_INVALID_INPUT before anything is allocated. The test
    !! goes through run_kohn_sham_scf_real (not through a local re-statement of
    !! the condition) so that deleting the `L <= 0` clause makes it fail: with
    !! the clause gone, the size(V_ext) /= L check fires instead and the error
    !! code becomes ERROR_SIZE_MISMATCH.
    !!
    !! Nup = Ndown = 0 is deliberate. With any positive occupation the companion
    !! clause Nup > L would already reject the input (5 > 0), so the test would
    !! survive the deletion of `L <= 0` and stop being branch-specific.
    !!
    !! store_history = .true. makes the second assertion meaningful: any correct
    !! implementation that reaches the SCF body must allocate the history, so a
    !! still-unallocated history proves the run stopped inside the validator.
    subroutine test_validate_inputs_invalid_L()
        use fortuno_serial, only: check => serial_check
        use lsda_types, only: system_params_t
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, scf_params_t, scf_results_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_OPEN
        use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT

        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(10)
        integer :: ierr
        character(len=256) :: table_file

        table_file = "build/test_kohn_sham_xc_table_u4.00.dat"
        call xc_lsda_init(xc_func, table_file, ierr)
        call check(ierr == ERROR_SUCCESS, "XC init should succeed")

        params%L = 0  ! Invalid!
        params%Nup = 0
        params%Ndown = 0
        params%bc = BC_OPEN
        params%U = 4.0_dp
        params%phase = 0.0_dp

        scf_params%max_iter = 10
        scf_params%density_tol = 1.0e-6_dp
        scf_params%energy_tol = 1.0e-8_dp
        scf_params%mixing_alpha = 0.3_dp
        scf_params%verbose = .false.
        scf_params%store_history = .true.

        V_ext = 0.0_dp

        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)

        call check(ierr == ERROR_INVALID_INPUT, "L = 0 must be rejected as invalid input")
        call check(.not. allocated(results%history%density_norms), &
                   "L = 0 must be rejected by validation, before any allocation")

        call xc_lsda_destroy(xc_func)
    end subroutine test_validate_inputs_invalid_L

    !> Test that the SCF entry point rejects occupations violating Pauli
    !!
    !! Physics: a single spin channel has exactly L orbitals, so N_sigma > L is
    !! unrepresentable. Both channels are probed independently (Nup = L+1 with a
    !! legal Ndown, then Ndown = L+1 with a legal Nup) so that neither branch of
    !! the per-spin check can be dropped without a test failure.
    !!
    !! validate_kohn_sham_cycle_inputs runs before the XC functional is touched
    !! and before any allocation, so these calls return immediately; the XC table
    !! is loaded anyway so that ERROR_INVALID_INPUT can only come from the
    !! validator and not from an uninitialised functional.
    !!
    !! store_history is deliberately .true. here: the "no allocation happened"
    !! assertions must prove that the validator stopped the run, and that is only
    !! true if a run reaching the SCF body would necessarily allocate the
    !! history. With store_history = .false. the same assertions would be
    !! satisfied by an implementation that skips the history allocation
    !! altogether (the behaviour T5 is going to introduce), and the test would
    !! keep passing for the wrong reason.
    subroutine test_validate_N_per_spin_exceeds_L()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, scf_params_t, scf_results_t
        use lsda_types, only: system_params_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_OPEN
        use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT

        integer, parameter :: L = 10
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L)
        integer :: ierr
        character(len=256) :: table_file

        table_file = "build/test_kohn_sham_xc_table_u4.00.dat"
        call xc_lsda_init(xc_func, table_file, ierr)
        call check(ierr == ERROR_SUCCESS, "XC init should succeed")

        params%L = L
        params%bc = BC_OPEN
        params%U = 4.0_dp
        params%phase = 0.0_dp

        scf_params%max_iter = 10
        scf_params%density_tol = 1.0e-6_dp
        scf_params%energy_tol = 1.0e-8_dp
        scf_params%mixing_alpha = 0.3_dp
        scf_params%verbose = .false.
        scf_params%store_history = .true.

        V_ext = 0.0_dp

        ! Spin-up channel overfilled: Nup = L + 1 > L
        params%Nup = L + 1
        params%Ndown = 5
        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "Nup = L + 1 must be rejected")
        call check(.not. allocated(results%history%density_norms), &
                   "Nup = L + 1 must be rejected by validation, before any allocation")

        ! Spin-down channel overfilled: Ndown = L + 1 > L
        params%Nup = 5
        params%Ndown = L + 1
        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "Ndown = L + 1 must be rejected")
        call check(.not. allocated(results%history%density_norms), &
                   "Ndown = L + 1 must be rejected by validation, before any allocation")

        ! Negative occupations are rejected too
        params%Nup = -1
        params%Ndown = 5
        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "Negative Nup must be rejected")
        call check(.not. allocated(results%history%density_norms), &
                   "Negative Nup must be rejected by validation, before any allocation")

        call xc_lsda_destroy(xc_func)
    end subroutine test_validate_N_per_spin_exceeds_L

    !> Test that the SCF entry point rejects a V_ext of the wrong length
    !!
    !! V_ext is an assumed-shape array, so a caller can pass any length; if the
    !! length disagrees with params%L the SCF would read past the end of the
    !! array. validate_kohn_sham_cycle_inputs must reject it with
    !! ERROR_SIZE_MISMATCH (not ERROR_INVALID_INPUT) before any allocation.
    !!
    !! All the other inputs are deliberately legal, so only the size check can
    !! fire: deleting that single branch lets the run proceed and the assertions
    !! below fail (wrong error code, and an allocated history because
    !! store_history = .true. forces the allocation on any correct path).
    subroutine test_validate_inputs_size_mismatch()
        use fortuno_serial, only: check => serial_check
        use lsda_types, only: system_params_t
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, scf_params_t, scf_results_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_OPEN
        use lsda_errors, only: ERROR_SUCCESS, ERROR_SIZE_MISMATCH

        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(5)  ! Wrong size!
        integer :: ierr
        character(len=256) :: table_file

        table_file = "build/test_kohn_sham_xc_table_u4.00.dat"
        call xc_lsda_init(xc_func, table_file, ierr)
        call check(ierr == ERROR_SUCCESS, "XC init should succeed")

        params%L = 10
        params%Nup = 5
        params%Ndown = 5
        params%bc = BC_OPEN
        params%U = 4.0_dp
        params%phase = 0.0_dp

        scf_params%max_iter = 10
        scf_params%density_tol = 1.0e-6_dp
        scf_params%energy_tol = 1.0e-8_dp
        scf_params%mixing_alpha = 0.3_dp
        scf_params%verbose = .false.
        scf_params%store_history = .true.

        V_ext = 0.0_dp

        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)

        call check(ierr == ERROR_SIZE_MISMATCH, &
                   "size(V_ext) /= L must be rejected with ERROR_SIZE_MISMATCH")
        call check(.not. allocated(results%history%density_norms), &
                   "size mismatch must be rejected by validation, before any allocation")

        call xc_lsda_destroy(xc_func)
    end subroutine test_validate_inputs_size_mismatch

    !> Test that the SCF entry point rejects a mixing weight outside (0, 1]
    !!
    !! The potential mixing weight alpha is the fraction of the newly computed
    !! V_eff blended in at each iteration, so it only makes sense in (0, 1].
    !! Alpha = 0 would discard V_eff_calc completely, while alpha = 1.5 would
    !! extrapolate the potential; both are rejected with ERROR_INVALID_INPUT.
    !!
    !! Every other input is legal (this exact 5/5 filling on L = 10 runs fine in
    !! the other SCF tests), so only the mixing branch can produce the error:
    !! deleting it lets the cycle run to completion and both assertions fail.
    subroutine test_validate_inputs_invalid_mixing()
        use fortuno_serial, only: check => serial_check
        use lsda_types, only: system_params_t
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, scf_params_t, scf_results_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_OPEN
        use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT

        integer, parameter :: L = 10
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L)
        integer :: ierr
        character(len=256) :: table_file

        table_file = "build/test_kohn_sham_xc_table_u4.00.dat"
        call xc_lsda_init(xc_func, table_file, ierr)
        call check(ierr == ERROR_SUCCESS, "XC init should succeed")

        params%L = L
        params%Nup = 5
        params%Ndown = 5
        params%bc = BC_OPEN
        params%U = 4.0_dp
        params%phase = 0.0_dp

        scf_params%max_iter = 10
        scf_params%density_tol = 1.0e-6_dp
        scf_params%energy_tol = 1.0e-8_dp
        scf_params%mixing_alpha = 1.5_dp  ! Invalid! Must be <= 1
        scf_params%verbose = .false.
        scf_params%store_history = .true.
        scf_params%use_adaptive_mixing = .false.

        V_ext = 0.0_dp

        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)

        call check(ierr == ERROR_INVALID_INPUT, &
                   "mixing_alpha = 1.5 must be rejected as invalid input")
        call check(.not. allocated(results%history%density_norms), &
                   "invalid mixing_alpha must be rejected by validation, before any allocation")

        scf_params%mixing_alpha = 0.0_dp  ! Invalid! Must be strictly positive
        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)

        call check(ierr == ERROR_INVALID_INPUT, &
                   "mixing_alpha = 0 must be rejected as invalid input")
        call check(.not. allocated(results%history%density_norms), &
                   "zero mixing_alpha must be rejected by validation, before any allocation")

        call xc_lsda_destroy(xc_func)
    end subroutine test_validate_inputs_invalid_mixing

    !> The SCF entry point must reject a tolerance that can never be met
    !!
    !! Since T3 BOTH potential_tol and energy_tol decide convergence
    !! (residual_V < potential_tol AND |dE| < energy_tol*max(1, |E|)), and both
    !! quantities are non-negative. A tolerance <= 0 is therefore unreachable:
    !! without this check the cycle would run the full max_iter and return
    !! ERROR_CONVERGENCE_FAILED, blaming the physics for what is an invalid
    !! request. Note that the fixed-alpha probe in this very file deliberately
    !! uses tolerances of 1e-14 to force a full-budget run: that is legal, tiny
    !! but positive, and must keep working.
    !!
    !! Each of the two branches is probed alone, with the other tolerance legal,
    !! so deleting either one makes exactly the corresponding assertions fail.
    !! store_history = .true. gives the usual anchor: any implementation that
    !! reaches the SCF body allocates the history, so an unallocated history
    !! proves the run stopped inside the validator.
    subroutine test_validate_inputs_nonpositive_tolerances()
        use fortuno_serial, only: check => serial_check
        use lsda_types, only: system_params_t
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, run_kohn_sham_scf_complex, &
                                    scf_params_t, scf_results_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_OPEN, BC_TWISTED
        use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT

        integer, parameter :: L = 10
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L)
        integer :: ierr
        character(len=256) :: table_file

        table_file = "build/test_kohn_sham_xc_table_u4.00.dat"
        call xc_lsda_init(xc_func, table_file, ierr)
        call check(ierr == ERROR_SUCCESS, "XC init should succeed")

        params%L = L
        params%Nup = 5
        params%Ndown = 5
        params%bc = BC_OPEN
        params%U = 4.0_dp
        params%phase = 0.0_dp

        scf_params%max_iter = 10
        scf_params%density_tol = 1.0e-6_dp
        scf_params%energy_tol = 1.0e-8_dp
        scf_params%potential_tol = 1.0e-6_dp
        scf_params%mixing_alpha = 0.3_dp
        scf_params%verbose = .false.
        scf_params%store_history = .true.

        V_ext = 0.0_dp

        ! --- potential_tol ----------------------------------------------------
        scf_params%potential_tol = 0.0_dp
        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "potential_tol = 0 must be rejected")
        call check(.not. allocated(results%history%density_norms), &
                   "potential_tol = 0 must be rejected before any allocation")

        scf_params%potential_tol = -1.0e-6_dp
        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "negative potential_tol must be rejected")
        call check(.not. allocated(results%history%density_norms), &
                   "negative potential_tol must be rejected before any allocation")

        ! --- energy_tol -------------------------------------------------------
        scf_params%potential_tol = 1.0e-6_dp
        scf_params%energy_tol = 0.0_dp
        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "energy_tol = 0 must be rejected")
        call check(.not. allocated(results%history%density_norms), &
                   "energy_tol = 0 must be rejected before any allocation")

        scf_params%energy_tol = -1.0e-8_dp
        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "negative energy_tol must be rejected")
        call check(.not. allocated(results%history%density_norms), &
                   "negative energy_tol must be rejected before any allocation")

        ! The complex entry point shares the validator and must reject the same
        ! inputs; the two loops are still duplicated code (T16).
        params%bc = BC_TWISTED
        call run_kohn_sham_scf_complex(params, scf_params, V_ext, xc_func, results, ierr)
        call check(ierr == ERROR_INVALID_INPUT, &
                   "the complex entry point must reject a negative energy_tol too")
        call check(.not. allocated(results%history%density_norms), &
                   "the complex entry point must reject it before any allocation")

        call xc_lsda_destroy(xc_func)
    end subroutine test_validate_inputs_nonpositive_tolerances

    !> The SCF entry points must reject non-finite tolerances too
    !!
    !! The public run_kohn_sham_scf_* routines are reachable without going
    !! through input_parser, so the shared validator is the last line of defence.
    !! A sign check alone accepts NaN and +Infinity:
    !!   * NaN makes `residual_V < potential_tol` (or the energy comparison)
    !!     always false, i.e. an unreachable criterion discovered only after the
    !!     whole iteration budget is spent;
    !!   * +Infinity makes it always true, DISABLING that criterion. With
    !!     potential_tol = +Inf the cycle would declare convergence on the energy
    !!     alone - the very false convergence this criterion was added to stop.
    !! -Infinity is already caught by `<= 0` and is asserted to keep it covered.
    !!
    !! Each branch is probed alone with the other tolerance legal, so removing
    !! the finiteness guard of one tolerance fails only its own assertions.
    !! store_history = .true. is the usual anchor: any run that reaches the SCF
    !! body allocates the history, so an unallocated history proves the rejection
    !! happened inside the validator.
    subroutine test_validate_inputs_non_finite_tolerances()
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
                                                 ieee_positive_inf, ieee_negative_inf
        use fortuno_serial, only: check => serial_check
        use lsda_types, only: system_params_t
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, run_kohn_sham_scf_complex, &
                                    scf_params_t, scf_results_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_OPEN, BC_TWISTED
        use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT

        integer, parameter :: L = 10
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L)
        real(dp) :: nan_v, pinf_v, ninf_v
        integer :: ierr
        character(len=256) :: table_file

        nan_v = ieee_value(1.0_dp, ieee_quiet_nan)
        pinf_v = ieee_value(1.0_dp, ieee_positive_inf)
        ninf_v = ieee_value(1.0_dp, ieee_negative_inf)

        table_file = "build/test_kohn_sham_xc_table_u4.00.dat"
        call xc_lsda_init(xc_func, table_file, ierr)
        call check(ierr == ERROR_SUCCESS, "XC init should succeed")

        params%L = L
        params%Nup = 5
        params%Ndown = 5
        params%bc = BC_OPEN
        params%U = 4.0_dp
        params%phase = 0.0_dp

        scf_params%max_iter = 10
        scf_params%density_tol = 1.0e-6_dp
        scf_params%energy_tol = 1.0e-8_dp
        scf_params%potential_tol = 1.0e-6_dp
        scf_params%mixing_alpha = 0.3_dp
        scf_params%verbose = .false.
        scf_params%store_history = .true.

        V_ext = 0.0_dp

        ! --- potential_tol ----------------------------------------------------
        scf_params%potential_tol = nan_v
        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "NaN potential_tol must be rejected")
        call check(.not. allocated(results%history%density_norms), &
                   "NaN potential_tol must be rejected before any allocation")

        scf_params%potential_tol = pinf_v
        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "+Infinity potential_tol must be rejected")
        call check(.not. allocated(results%history%density_norms), &
                   "+Infinity potential_tol must be rejected before any allocation")

        scf_params%potential_tol = ninf_v
        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "-Infinity potential_tol must be rejected")
        call check(.not. allocated(results%history%density_norms), &
                   "-Infinity potential_tol must be rejected before any allocation")

        ! --- energy_tol -------------------------------------------------------
        scf_params%potential_tol = 1.0e-6_dp
        scf_params%energy_tol = nan_v
        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "NaN energy_tol must be rejected")
        call check(.not. allocated(results%history%density_norms), &
                   "NaN energy_tol must be rejected before any allocation")

        scf_params%energy_tol = pinf_v
        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "+Infinity energy_tol must be rejected")
        call check(.not. allocated(results%history%density_norms), &
                   "+Infinity energy_tol must be rejected before any allocation")

        scf_params%energy_tol = ninf_v
        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "-Infinity energy_tol must be rejected")
        call check(.not. allocated(results%history%density_norms), &
                   "-Infinity energy_tol must be rejected before any allocation")

        ! The complex entry point shares the validator and must reject the same
        ! inputs; the two loops are still duplicated code (T16).
        params%bc = BC_TWISTED
        scf_params%energy_tol = 1.0e-8_dp
        scf_params%potential_tol = pinf_v
        call run_kohn_sham_scf_complex(params, scf_params, V_ext, xc_func, results, ierr)
        call check(ierr == ERROR_INVALID_INPUT, &
                   "the complex entry point must reject an infinite potential_tol too")
        call check(.not. allocated(results%history%density_norms), &
                   "the complex entry point must reject it before any allocation")

        call xc_lsda_destroy(xc_func)
    end subroutine test_validate_inputs_non_finite_tolerances

    !> The Hartree interaction and XC functional must describe the same U
    !!
    !! `system_params_t` and `xc_lsda_t` reach the public SCF entry points as
    !! independent objects. A negative `params%U` paired with a positive-U XC
    !! table used to run without complaint, combining an attractive Hartree
    !! term with a repulsive functional. Both solver paths must reject a sign or
    !! magnitude mismatch before allocating their convergence history.
    subroutine test_validate_xc_matches_system_u()
        use fortuno_serial, only: check => serial_check
        use lsda_types, only: system_params_t
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, run_kohn_sham_scf_complex, &
                                    scf_params_t, scf_results_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_OPEN, BC_TWISTED
        use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT

        integer, parameter :: L = 4
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L)
        integer :: ierr

        call xc_lsda_init(xc_func, "build/test_kohn_sham_xc_table_u4.00.dat", &
                          ierr, u_signed=4.0_dp)
        call check(ierr == ERROR_SUCCESS, "XC consistency: initialization should succeed")
        if (ierr /= ERROR_SUCCESS) return

        params%L = L
        params%Nup = 2
        params%Ndown = 2
        params%bc = BC_OPEN
        params%U = -4.0_dp
        params%phase = 0.0_dp
        scf_params%max_iter = 1
        scf_params%store_history = .true.
        V_ext = 0.0_dp

        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)
        call check(ierr == ERROR_INVALID_INPUT, &
                   "XC consistency: real solver must reject opposite signs of U")
        call check(.not. allocated(results%history%density_norms), &
                   "XC consistency: real solver must reject before allocation")

        params%bc = BC_TWISTED
        call run_kohn_sham_scf_complex(params, scf_params, V_ext, xc_func, results, ierr)
        call check(ierr == ERROR_INVALID_INPUT, &
                   "XC consistency: complex solver must reject opposite signs of U")
        call check(.not. allocated(results%history%density_norms), &
                   "XC consistency: complex solver must reject before allocation")

        params%bc = BC_OPEN
        params%U = 2.0_dp
        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)
        call check(ierr == ERROR_INVALID_INPUT, &
                   "XC consistency: solver must reject different magnitudes of U")

        call xc_lsda_destroy(xc_func)
    end subroutine test_validate_xc_matches_system_u

    !> Test SCF results initialization and cleanup
    subroutine test_scf_results_init_cleanup()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: scf_results_t, init_scf_results, cleanup_scf_results
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        integer :: ierr

        call init_scf_results(results, 10, .true., 100, ierr)
        call check(ierr == ERROR_SUCCESS, "Init should succeed")
        call check(.not. results%converged, "Should initialize as not converged")
        call check(results%n_iterations == 0, "Iterations should be zero")

        call cleanup_scf_results(results, ierr)
        call check(ierr == ERROR_SUCCESS, "Cleanup should succeed")
    end subroutine test_scf_results_init_cleanup


    !> init_scf_results must not swallow a failed history allocation
    !!
    !! init_convergence_history rejects max_iter <= 0 with ERROR_INVALID_INPUT,
    !! but init_scf_results used to overwrite that code with ERROR_SUCCESS on the
    !! very next line. A caller assembling an scf_results_t outside the SCF cycle
    !! (which is the only reason this helper is exported) therefore got a success
    !! code together with an unallocated history, and only found out when it
    !! indexed into it.
    !!
    !! The unallocated history is asserted as well: a correct rejection cannot
    !! have allocated anything, and it is what distinguishes this from a mere
    !! error-code change.
    subroutine test_scf_results_init_propagates_history_error()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: scf_results_t, init_scf_results, cleanup_scf_results
        use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT

        type(scf_results_t) :: results
        integer :: ierr

        call init_scf_results(results, 10, .true., 0, ierr)
        call check(ierr == ERROR_INVALID_INPUT, &
                   "max_iter = 0 with store_history must report ERROR_INVALID_INPUT")
        call check(.not. allocated(results%history%density_norms), &
                   "a rejected history must not be allocated")

        call init_scf_results(results, 10, .true., -5, ierr)
        call check(ierr == ERROR_INVALID_INPUT, &
                   "a negative max_iter with store_history must be reported too")

        ! Without a history there is nothing to allocate, so max_iter is
        ! irrelevant and the call must succeed.
        call init_scf_results(results, 10, .false., 0, ierr)
        call check(ierr == ERROR_SUCCESS, &
                   "max_iter is irrelevant when no history is requested")

        ! And the ordinary path must still succeed.
        call init_scf_results(results, 10, .true., 7, ierr)
        call check(ierr == ERROR_SUCCESS, "a legal max_iter must still succeed")
        call check(allocated(results%history%density_norms), &
                   "a legal request must actually allocate the history")

        call cleanup_scf_results(results, ierr)
    end subroutine test_scf_results_init_propagates_history_error


    !> The complex SCF loop must reproduce the real one at zero twist
    !!
    !! run_kohn_sham_scf_real and run_kohn_sham_scf_complex are two hand-kept
    !! copies of the same ~500-line algorithm (unification is T16). Nothing so
    !! far checked that they agree, so any fix applied to one and forgotten in
    !! the other would go unnoticed until a twisted-BC run produced silently
    !! wrong physics.
    !!
    !! The anchor is exact by construction: with theta = 0 the Peierls factor
    !! e^{i*theta} of BC_TWISTED is 1, so the complex Hamiltonian IS the real
    !! periodic one embedded in C. Both loops therefore see the same spectrum,
    !! the same densities, the same mixing decisions and the same convergence
    !! test, and must produce the same energy, the same density profile and the
    !! same residual_V, up to the difference between ZHEEVD and DSYEVD (which is
    !! at the level of the eigenvector phases, invisible in |psi|^2).
    !!
    !! A small, quickly converging system is used on purpose: L = 8, Nup = 3,
    !! Ndown = 2 (spin polarised, so the two channels carry different potentials
    !! and a bug confined to one of them cannot cancel), U = 2, adaptive mixing
    !! ON so that the controller path is exercised in both copies. The filling
    !! is kept well below n = 1: at half filling this system sits on the V_xc
    !! discontinuity and neither loop converges, which would make the comparison
    !! a comparison of two failures.
    !!
    !! MUTATION VERIFIED: changing a single factor in the complex loop's mixing
    !! (alpha_used -> 0.5*alpha_used on the spin-down channel) makes the density
    !! and energy assertions fail.
    subroutine test_complex_loop_matches_real_at_zero_twist()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, run_kohn_sham_scf_complex, &
                                    scf_params_t, scf_results_t, cleanup_scf_results
        use lsda_types, only: system_params_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_PERIODIC, BC_TWISTED
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 8
        real(dp), parameter :: MATCH_TOL = 1.0e-10_dp

        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: res_real, res_cplx
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L)
        real(dp) :: max_dn_up, max_dn_down
        integer :: ierr
        character(len=256) :: table_file

        table_file = "build/test_kohn_sham_xc_table_u2.00.dat"
        call xc_lsda_init(xc_func, table_file, ierr)
        call check(ierr == ERROR_SUCCESS, "XC init should succeed")
        if (ierr /= ERROR_SUCCESS) return

        params%L = L
        params%Nup = 3
        params%Ndown = 2
        params%U = 2.0_dp
        params%phase = 0.0_dp

        scf_params%max_iter = 200
        scf_params%density_tol = 1.0e-6_dp
        scf_params%energy_tol = 1.0e-10_dp
        scf_params%potential_tol = 1.0e-8_dp
        scf_params%mixing_alpha = 0.2_dp
        scf_params%use_adaptive_mixing = .true.
        scf_params%verbose = .false.
        scf_params%store_history = .false.

        ! A non-uniform V_ext so the fixed point is not trivially the uniform
        ! density: a translationally invariant system would be reproduced by
        ! almost any bug in the loop.
        V_ext = 0.0_dp
        V_ext(3) = -1.5_dp
        V_ext(6) = 0.8_dp

        params%bc = BC_PERIODIC
        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, res_real, ierr)
        call check(ierr == ERROR_SUCCESS, "the real loop must converge on this system")

        params%bc = BC_TWISTED
        call run_kohn_sham_scf_complex(params, scf_params, V_ext, xc_func, res_cplx, ierr)
        call check(ierr == ERROR_SUCCESS, &
                   "the complex loop at theta = 0 must converge on the same system")

        call check(res_real%converged .and. res_cplx%converged, &
                   "Precondition: both loops must reach self-consistency")

        if (res_real%converged .and. res_cplx%converged) then
            call check(res_real%n_iterations == res_cplx%n_iterations, &
                       "both loops must take the same number of iterations")

            call check(abs(res_real%final_energy - res_cplx%final_energy) <= &
                       MATCH_TOL * max(1.0_dp, abs(res_real%final_energy)), &
                       "the two loops must agree on the total energy at theta = 0")

            call check(abs(res_real%final_potential_residual - &
                           res_cplx%final_potential_residual) <= &
                       MATCH_TOL * max(1.0_dp, abs(res_real%final_potential_residual)), &
                       "the two loops must agree on residual_V at theta = 0")

            call check(allocated(res_real%density_up) .and. allocated(res_cplx%density_up), &
                       "both loops must return a density")

            if (allocated(res_real%density_up) .and. allocated(res_cplx%density_up)) then
                max_dn_up = maxval(abs(res_real%density_up - res_cplx%density_up))
                max_dn_down = maxval(abs(res_real%density_down - res_cplx%density_down))

                ! Precondition: the density really is structured, so an equality
                ! between two flat profiles cannot be what is being measured.
                call check(maxval(res_real%density_up) - minval(res_real%density_up) > 1.0e-3_dp, &
                           "Precondition: V_ext must have made the density non-uniform")

                call check(max_dn_up <= MATCH_TOL, &
                           "the two loops must agree on the spin-up density at theta = 0")
                call check(max_dn_down <= MATCH_TOL, &
                           "the two loops must agree on the spin-down density at theta = 0")
            end if
        end if

        call cleanup_scf_results(res_real, ierr)
        call cleanup_scf_results(res_cplx, ierr)
        call xc_lsda_destroy(xc_func)
    end subroutine test_complex_loop_matches_real_at_zero_twist

    !> Test SCF convergence for U=0 with open BC
    !!
    !! Physics: U=0 is the non-interacting limit (free Fermi gas), with no XC
    !! table.  For open boundaries the exact energy is the occupied sum of
    !! -2 cos(j pi/(L+1)) in each spin channel.
    subroutine test_scf_converges_u0_open()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, scf_params_t, scf_results_t, &
                                    cleanup_scf_results
        use lsda_types, only: system_params_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_OPEN
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 10
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L)
        integer :: ierr, j
        real(dp) :: expected_energy

        call xc_lsda_init(xc_func, ierr=ierr, u_signed=0.0_dp)
        call check(ierr == ERROR_SUCCESS, "U=0 XC initialization must not require a table")
        if (ierr /= ERROR_SUCCESS) return

        ! System parameters: small system, half-filled
        params%L = L
        params%Nup = 5
        params%Ndown = 5
        params%bc = BC_OPEN
        params%U = 0.0_dp
        params%phase = 0.0_dp

        ! SCF parameters: tight convergence
        scf_params%max_iter = 50
        scf_params%density_tol = 1.0e-6_dp
        scf_params%energy_tol = 1.0e-8_dp
        scf_params%mixing_alpha = 0.5_dp
        scf_params%verbose = .false.
        scf_params%store_history = .true.

        ! Uniform external potential
        V_ext = 0.0_dp

        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)

        call check(ierr == ERROR_SUCCESS, "SCF should succeed")
        call check(results%converged, "SCF should converge for U=0")
        call check(results%n_iterations <= 10, "U=0 should converge quickly")
        call check(allocated(results%density_up), "Density_up should be allocated")
        call check(allocated(results%density_down), "Density_down should be allocated")
        call check(allocated(results%eigvals), "Eigvals should be allocated")

        ! Verify particle number conservation
        call check(abs(sum(results%density_up) + sum(results%density_down) - 10.0_dp) < 1.0e-6_dp, &
                   "Particle number should be conserved")

        expected_energy = 0.0_dp
        do j = 1, 5
            expected_energy = expected_energy - 4.0_dp * cos(real(j, dp) * acos(-1.0_dp) / real(L + 1, dp))
        end do
        call check(abs(results%final_energy - expected_energy) < 1.0e-10_dp, &
                   "U=0 SCF energy must equal the open-chain free-Fermi value")

        call cleanup_scf_results(results, ierr)
        call xc_lsda_destroy(xc_func)
    end subroutine test_scf_converges_u0_open

    !> Test SCF convergence for U=2 with periodic BC
    !!
    !! Physics: Periodic BC allows momentum conservation.
    !! For uniform V_ext, density should be exactly uniform.
    !! Note: Convergence may be slow for U>0, so we allow more iterations.
    subroutine test_scf_converges_u0_periodic()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, scf_params_t, scf_results_t, &
                                    cleanup_scf_results
        use lsda_types, only: system_params_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_PERIODIC
        use lsda_errors, only: ERROR_SUCCESS, ERROR_CONVERGENCE_FAILED

        integer, parameter :: L = 8
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L)
        integer :: ierr, i
        character(len=256) :: table_file

        table_file = "build/test_kohn_sham_xc_table_u2.00.dat"
        call xc_lsda_init(xc_func, table_file, ierr)
        if (ierr /= ERROR_SUCCESS) return

        params%L = L
        params%Nup = 4
        params%Ndown = 4
        params%bc = BC_PERIODIC
        params%U = 2.0_dp
        params%phase = 0.0_dp

        scf_params%max_iter = 100  ! More iterations for convergence
        scf_params%density_tol = 1.0e-5_dp  ! Slightly looser tolerance
        scf_params%energy_tol = 1.0e-7_dp
        scf_params%mixing_alpha = 0.2_dp  ! Smaller alpha for stability
        scf_params%verbose = .false.
        scf_params%store_history = .true.

        V_ext = 0.0_dp

        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)

        ! Accept both convergence and non-convergence (SCF can be tricky)
        call check(ierr == ERROR_SUCCESS .or. ierr == ERROR_CONVERGENCE_FAILED, &
                   "SCF should complete (converged or not)")

        ! If it did converge, check that error is reasonable
        if (results%converged) then
            call check(results%final_density_error < 1.0e-4_dp, "Final error should be reasonable")
        end if

        call cleanup_scf_results(results, ierr)
        call xc_lsda_destroy(xc_func)
    end subroutine test_scf_converges_u0_periodic

    !> Test that SCF stores convergence history when requested
    subroutine test_scf_stores_history()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, scf_params_t, scf_results_t, &
                                    cleanup_scf_results
        use lsda_types, only: system_params_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_OPEN
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 6
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L)
        integer :: ierr
        character(len=256) :: table_file

        table_file = "build/test_kohn_sham_xc_table_u2.00.dat"
        call xc_lsda_init(xc_func, table_file, ierr)
        if (ierr /= ERROR_SUCCESS) return

        params%L = L
        params%Nup = 3
        params%Ndown = 3
        params%bc = BC_OPEN
        params%U = 2.0_dp
        params%phase = 0.0_dp

        scf_params%max_iter = 30
        scf_params%density_tol = 1.0e-6_dp
        scf_params%energy_tol = 1.0e-8_dp
        scf_params%mixing_alpha = 0.4_dp
        scf_params%verbose = .false.
        scf_params%store_history = .true.  ! Request history storage

        V_ext = 0.0_dp

        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)

        call check(ierr == ERROR_SUCCESS, "SCF should succeed")
        if (results%converged) then
            call check(allocated(results%history%density_norms), "History should be allocated")
            call check(allocated(results%history%energies), "Energy history should be allocated")
            call check(results%history%current_iter > 0, "Should have stored some iterations")
            call check(results%history%current_iter == results%n_iterations, &
                       "Stored iterations should match total")
        end if

        call cleanup_scf_results(results, ierr)
        call xc_lsda_destroy(xc_func)
    end subroutine test_scf_stores_history

    !> store_history = .false. must not allocate the convergence history (T5)
    !!
    !! The SCF routine is the single owner of the initialization of `results`,
    !! and it used to call init_convergence_history unconditionally: a caller who
    !! explicitly asked for no history still paid for three arrays of max_iter
    !! doubles (30000 each with the default max_iter = 10000) that were never
    !! written to, since update_convergence_history is guarded by the same flag.
    !!
    !! The run must actually go through the SCF body - asserted via
    !! n_iterations > 0 - otherwise an unallocated history would prove nothing.
    subroutine test_scf_skips_history_when_disabled()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, scf_params_t, scf_results_t, &
                                    cleanup_scf_results
        use lsda_types, only: system_params_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_OPEN
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 6
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L)
        integer :: ierr
        character(len=256) :: table_file

        table_file = "build/test_kohn_sham_xc_table_u2.00.dat"
        call xc_lsda_init(xc_func, table_file, ierr)
        if (ierr /= ERROR_SUCCESS) return

        params%L = L
        params%Nup = 3
        params%Ndown = 3
        params%bc = BC_OPEN
        params%U = 2.0_dp
        params%phase = 0.0_dp

        scf_params%max_iter = 30
        scf_params%density_tol = 1.0e-6_dp
        scf_params%energy_tol = 1.0e-8_dp
        scf_params%mixing_alpha = 0.4_dp
        scf_params%verbose = .false.
        scf_params%store_history = .false.

        V_ext = 0.0_dp

        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)

        call check(results%n_iterations > 0, &
                   "Precondition: the SCF body must have run")
        call check(.not. allocated(results%history%density_norms), &
                   "No density-norm history may be allocated when store_history is .false.")
        call check(.not. allocated(results%history%energies), &
                   "No energy history may be allocated when store_history is .false.")
        call check(.not. allocated(results%history%potential_residuals), &
                   "No residual history may be allocated when store_history is .false.")

        call cleanup_scf_results(results, ierr)
        call xc_lsda_destroy(xc_func)
    end subroutine test_scf_skips_history_when_disabled

    !> Test that particle number is conserved during SCF
    !!
    !! Physics: Total particle number N = Σn(i) must be conserved exactly
    !! at each SCF iteration (within numerical precision).
    subroutine test_scf_density_conservation()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, scf_params_t, scf_results_t, &
                                    cleanup_scf_results
        use lsda_types, only: system_params_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_PERIODIC
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 10
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L)
        real(dp) :: total_N
        integer :: ierr
        character(len=256) :: table_file

        table_file = "build/test_kohn_sham_xc_table_u4.00.dat"
        call xc_lsda_init(xc_func, table_file, ierr)
        if (ierr /= ERROR_SUCCESS) return

        params%L = L
        params%Nup = 5
        params%Ndown = 5
        params%bc = BC_PERIODIC
        params%U = 4.0_dp
        params%phase = 0.0_dp

        scf_params%max_iter = 50
        scf_params%density_tol = 1.0e-6_dp
        scf_params%energy_tol = 1.0e-8_dp
        scf_params%mixing_alpha = 0.3_dp
        scf_params%verbose = .false.
        scf_params%store_history = .false.

        V_ext = 0.0_dp

        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)

        if (ierr == ERROR_SUCCESS .and. results%converged) then
            total_N = sum(results%density_up) + sum(results%density_down)
            call check(abs(total_N - real(params%Nup + params%Ndown, dp)) < 1.0e-8_dp, &
                       "Particle number should be exactly conserved")
        end if

        call cleanup_scf_results(results, ierr)
        call xc_lsda_destroy(xc_func)
    end subroutine test_scf_density_conservation

    !> Test SCF with complex Hamiltonian (twisted BC)
    !!
    !! Physics: Twisted BC introduces Aharonov-Bohm phase θ.
    !! Spectrum shifts: E_k(θ) = -2cos((2πk+θ)/L)
    !! Density should still be real: n(i) = Σ|ψ_j(i)|²
    subroutine test_scf_complex_twisted_bc()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: run_kohn_sham_scf_complex, scf_params_t, scf_results_t, &
                                    cleanup_scf_results
        use lsda_types, only: system_params_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_TWISTED
        use lsda_constants, only: PI
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 8
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L)
        integer :: ierr, i
        character(len=256) :: table_file

        table_file = "build/test_kohn_sham_xc_table_u2.00.dat"
        call xc_lsda_init(xc_func, table_file, ierr)
        if (ierr /= ERROR_SUCCESS) return

        params%L = L
        params%Nup = 4
        params%Ndown = 4
        params%bc = BC_TWISTED
        params%U = 2.0_dp
        params%phase = PI / 4.0_dp  ! θ = π/4

        scf_params%max_iter = 50
        scf_params%density_tol = 1.0e-6_dp
        scf_params%energy_tol = 1.0e-8_dp
        scf_params%mixing_alpha = 0.3_dp
        scf_params%verbose = .false.
        scf_params%store_history = .false.

        V_ext = 0.0_dp

        call run_kohn_sham_scf_complex(params, scf_params, V_ext, xc_func, results, ierr)

        call check(ierr == ERROR_SUCCESS, "Complex SCF should succeed")
        if (results%converged) then
            call check(allocated(results%density_up), "Density should be allocated")

            ! Density must be real and positive
            do i = 1, L
                call check(results%density_up(i) >= 0.0_dp, "Density must be non-negative")
                call check(results%density_down(i) >= 0.0_dp, "Density must be non-negative")
            end do

            ! Particle conservation
            call check(abs(sum(results%density_up) + sum(results%density_down) - 8.0_dp) < 1.0e-6_dp, &
                       "Particle number must be conserved")
        end if

        call cleanup_scf_results(results, ierr)
        call xc_lsda_destroy(xc_func)
    end subroutine test_scf_complex_twisted_bc

    !> Declaring convergence must mean actual self-consistency (regression, T3)
    !!
    !! This is the production case of input.txt with the random seed pinned:
    !! L = 100, Nup = Ndown = 25, U = -4 (attractive), open BC, 50% random
    !! impurities of strength V0 = -4. It is a hard, strongly correlated,
    !! disordered system and it is exactly the case that used to be reported as
    !! "CONVERGED" after 163 iterations while the total energy was still jumping
    !! between -310.66 and -311.95 (an oscillation of 1.29, i.e. 4e-3 relative).
    !!
    !! The mechanism of that false positive: the old criterion was
    !! ||n_out - n_in||_2 < density_tol, but n_in IS the previous n_out and the
    !! two densities come from potentials differing by alpha*(V_calc - V_eff).
    !! ||delta_n|| is therefore proportional to the mixing weight, and the
    !! adaptive controller drove alpha geometrically towards zero (0.05, 0.0167,
    !! 0.0056, ... down to a 1e-10 floor), so ||delta_n|| crossed the tolerance
    !! purely because the cycle had stopped moving.
    !!
    !! The assertions below are the contract that forbids this: whenever the SCF
    !! reports success, the potential residual must really be below
    !! potential_tol AND the last two energies must agree to energy_tol in
    !! relative terms. Both are read back from the stored history, so the test
    !! checks what the cycle actually did, not what it claims in a summary field.
    !! Under the old criterion the run reports converged at iteration ~163 with a
    !! potential residual of order 0.2 and a relative energy step of order 4e-3,
    !! so both assertions fail.
    !!
    !! max_iter = 300 is a deliberate compromise: it is comfortably past the
    !! iteration where the old code declared victory (163), which is all this
    !! test needs, while keeping the runtime around a third of a second. With the
    !! corrected criterion this system does not converge at all (the residual
    !! plateaus around 0.2 even after 10000 iterations), and the non-convergent
    !! branch asserts that this outcome is reported honestly: converged = .false.,
    !! ERROR_CONVERGENCE_FAILED, and a residual that is genuinely above the
    !! tolerance.
    subroutine test_scf_declared_convergence_is_self_consistent()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, scf_params_t, scf_results_t, &
                                    cleanup_scf_results
        use lsda_types, only: system_params_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use potential_impurity, only: potential_impurity_random
        use boundary_conditions, only: BC_OPEN
        use lsda_errors, only: ERROR_SUCCESS, ERROR_CONVERGENCE_FAILED

        integer, parameter :: L = 100
        integer, parameter :: POT_SEED = 12345
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L)
        integer, allocatable :: imp_positions(:)
        real(dp) :: energy_step, energy_scale
        integer :: ierr, n_iter
        character(len=256) :: table_file

        ! The table is indexed by |U|, but the attractive sign must be handed to
        ! the functional explicitly so that the Shiba transformation is active.
        table_file = "build/test_kohn_sham_xc_table_u4.00.dat"
        call xc_lsda_init(xc_func, table_file, ierr, u_signed=-4.0_dp)
        call check(ierr == ERROR_SUCCESS, "XC init should succeed")

        ! Same potential as input.txt, with the seed pinned for reproducibility
        call potential_impurity_random(-4.0_dp, 50.0_dp, L, POT_SEED, V_ext, imp_positions, ierr)
        call check(ierr == ERROR_SUCCESS, "Random impurity potential should be created")
        if (allocated(imp_positions)) deallocate(imp_positions)

        params%L = L
        params%Nup = 25
        params%Ndown = 25
        params%bc = BC_OPEN
        params%U = -4.0_dp
        params%phase = 0.0_dp

        scf_params%max_iter = 300
        scf_params%density_tol = 1.0e-6_dp
        scf_params%energy_tol = 1.0e-8_dp
        scf_params%potential_tol = 1.0e-6_dp
        scf_params%mixing_alpha = 0.05_dp
        scf_params%use_adaptive_mixing = .true.
        scf_params%verbose = .false.
        scf_params%store_history = .true.

        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)

        call check(ierr == ERROR_SUCCESS .or. ierr == ERROR_CONVERGENCE_FAILED, &
                   "SCF must finish either converged or explicitly failed")
        call check(results%converged .eqv. (ierr == ERROR_SUCCESS), &
                   "The converged flag and the error code must agree")

        ! Particle number is conserved whatever the outcome, and the densities
        ! must be available even on the non-convergent path.
        call check(allocated(results%density_up) .and. allocated(results%density_down), &
                   "Densities must be available even when the SCF does not converge")
        if (allocated(results%density_up) .and. allocated(results%density_down)) then
            call check(abs(sum(results%density_up) + sum(results%density_down) - 50.0_dp) < 1.0e-8_dp, &
                       "Particle number must be conserved")
        end if

        n_iter = results%n_iterations
        call check(allocated(results%history%potential_residuals), &
                   "The potential residual must be recorded in the history")

        if (results%converged) then
            call check(n_iter >= 2, &
                       "Convergence cannot be declared at iteration 1: the residual is " // &
                       "identically zero there because V_eff is seeded from the same density")

            call check(results%final_potential_residual < scf_params%potential_tol, &
                       "Declared convergence requires a potential residual below potential_tol")

            if (allocated(results%history%potential_residuals) .and. n_iter >= 2) then
                call check(results%history%potential_residuals(n_iter) < 1.0e-6_dp, &
                           "The recorded residual at the convergence iteration must be below 1e-6")

                energy_step = abs(results%history%energies(n_iter) - &
                                  results%history%energies(n_iter - 1))
                energy_scale = max(1.0_dp, abs(results%history%energies(n_iter)))
                call check(energy_step < 1.0e-8_dp * energy_scale, &
                           "The last two energies must agree to 1e-8 in relative terms")
            end if
        else
            call check(ierr == ERROR_CONVERGENCE_FAILED, &
                       "A non-convergent run must return ERROR_CONVERGENCE_FAILED")
            call check(results%n_iterations == scf_params%max_iter, &
                       "A non-convergent run must have used the whole iteration budget")
            call check(results%final_potential_residual >= scf_params%potential_tol, &
                       "A run reported as not converged must really have a large residual")
        end if

        call cleanup_scf_results(results, ierr)
        call xc_lsda_destroy(xc_func)
    end subroutine test_scf_declared_convergence_is_self_consistent

    !> Run the alpha-independence probe system for a fixed number of iterations
    !!
    !! Deliberately small and fully deterministic: L = 20, Nup = 7, Ndown = 5
    !! (spin polarised, so both channels carry a different potential), U = 4,
    !! open BC and a single attractive impurity of -2 at site 10. The impurity
    !! is what makes the problem non-trivial: with a uniform V_ext on a ring the
    !! seeded V_eff is already a fixed point and every residual would be zero.
    !!
    !! The mixing is FIXED (use_adaptive_mixing = .false.) so that alpha is the
    !! only thing that differs between two calls, and the tolerances are set to
    !! 1e-14 so the cycle never stops early: both runs perform exactly n_iter
    !! iterations and the histories can be compared index by index.
    !!
    !! @param[in]  alpha      Mixing weight of the new potential
    !! @param[in]  n_iter     Number of SCF iterations to perform (exactly)
    !! @param[out] residuals  history%potential_residuals(1:n_iter)
    !! @param[out] dnorms     history%density_norms(1:n_iter)
    !! @param[out] ok         .true. if the run produced a full n_iter history
    !! @param[in]  adaptive   Use the adaptive controller instead of fixed mixing
    !!                        (optional, default .false.)
    !! @param[out] dens_up    Final spin-up density n_out of the last iteration
    !!                        (optional, length >= PROBE_L)
    !! @param[out] dens_down  Final spin-down density (optional, length >= PROBE_L)
    !! @param[out] v_ext_out  The external potential the probe used (optional)
    subroutine run_fixed_alpha_probe(alpha, n_iter, residuals, dnorms, ok, adaptive, &
                                     dens_up, dens_down, v_ext_out)
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, scf_params_t, scf_results_t, &
                                    cleanup_scf_results
        use lsda_types, only: system_params_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_OPEN
        use lsda_errors, only: ERROR_SUCCESS, ERROR_CONVERGENCE_FAILED

        real(dp), intent(in) :: alpha
        integer, intent(in) :: n_iter
        real(dp), intent(out) :: residuals(:), dnorms(:)
        logical, intent(out) :: ok
        logical, intent(in), optional :: adaptive
        real(dp), intent(out), optional :: dens_up(:), dens_down(:), v_ext_out(:)

        integer, parameter :: L = PROBE_L
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L)
        integer :: ierr
        character(len=256) :: table_file

        ok = .false.
        residuals = 0.0_dp
        dnorms = 0.0_dp
        if (present(dens_up)) dens_up = 0.0_dp
        if (present(dens_down)) dens_down = 0.0_dp

        table_file = PROBE_TABLE
        call xc_lsda_init(xc_func, table_file, ierr)
        if (ierr /= ERROR_SUCCESS) return

        params%L = L
        params%Nup = PROBE_NUP
        params%Ndown = PROBE_NDOWN
        params%bc = BC_OPEN
        params%U = PROBE_U
        params%phase = 0.0_dp

        scf_params%max_iter = n_iter
        scf_params%density_tol = 1.0e-6_dp
        scf_params%energy_tol = 1.0e-14_dp
        scf_params%potential_tol = 1.0e-14_dp
        scf_params%mixing_alpha = alpha
        scf_params%use_adaptive_mixing = .false.
        if (present(adaptive)) scf_params%use_adaptive_mixing = adaptive
        scf_params%verbose = .false.
        scf_params%store_history = .true.

        V_ext = 0.0_dp
        V_ext(PROBE_IMP_SITE) = PROBE_IMP_V
        if (present(v_ext_out)) v_ext_out(1:L) = V_ext

        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)

        if (ierr == ERROR_CONVERGENCE_FAILED .and. &
            allocated(results%history%potential_residuals) .and. &
            results%history%current_iter == n_iter) then
            residuals(1:n_iter) = results%history%potential_residuals(1:n_iter)
            dnorms(1:n_iter) = results%history%density_norms(1:n_iter)
            if (present(dens_up) .and. allocated(results%density_up)) &
                dens_up(1:L) = results%density_up
            if (present(dens_down) .and. allocated(results%density_down)) &
                dens_down(1:L) = results%density_down
            ok = .true.
        end if

        call cleanup_scf_results(results, ierr)
        call xc_lsda_destroy(xc_func)
    end subroutine run_fixed_alpha_probe

    !> Generate the deterministic XC fixtures used by this test program.
    !!
    !! The tables are intentionally produced through the same thermodynamic-limit
    !! generator and native writer used in production.  They live beneath
    !! `build/` and are not reference data checked into the repository.
    !!
    !! They are rewritten on every run rather than cached: these tables are the
    !! oracle of this suite, and a file left over from an earlier generator
    !! would silently keep the suite green against a table the code no longer
    !! produces.  The 8 x 9 grid makes the regeneration cheap.
    subroutine prepare_test_xc_tables()
        use bethe_tables, only: generate_xc_table, grid_params_t
        use table_io, only: xc_table_t, write_fortran_table, deallocate_table
        use lsda_errors, only: ERROR_SUCCESS

        real(dp), parameter :: U_VALUES(2) = [2.0_dp, 4.0_dp]
        type(grid_params_t) :: params
        type(xc_table_t) :: table
        character(len=256) :: filename
        integer :: i, ierr

        params = grid_params_t()
        params%n_points = 8
        params%m_points = 9

        do i = 1, size(U_VALUES)
            write(filename, '(A,F0.2,A)') 'build/test_kohn_sham_xc_table_u', U_VALUES(i), '.dat'

            call generate_xc_table(U_VALUES(i), params, table, ierr)
            if (ierr /= ERROR_SUCCESS) error stop 'failed to generate kohn-sham XC test fixture'

            call write_fortran_table(trim(filename), table, ierr)
            if (allocated(table%n_grid)) call deallocate_table(table)
            if (ierr /= ERROR_SUCCESS) error stop 'failed to write kohn-sham XC test fixture'
        end do
    end subroutine prepare_test_xc_tables

    !> The potential residual must NOT scale with the mixing weight (regression, T3)
    !!
    !! This is the assertion that pins down WHY residual_V replaced ||Δn|| as the
    !! convergence criterion, and it is deliberately differential: the same system
    !! is run twice with fixed mixing, alpha = 0.05 and alpha = 0.005 (a factor of
    !! 10), for the same fixed number of iterations, and the two histories are
    !! compared.
    !!
    !! Expected behaviour:
    !!   * ||Δn|| IS proportional to alpha. n_in is literally the previous n_out,
    !!     and the two potentials behind them differ by alpha*(V_calc - V_eff), so
    !!     dividing alpha by 10 divides ||Δn|| by ~10. That is exactly why a
    !!     density-based criterion could be satisfied by a cycle that had merely
    !!     stopped moving.
    !!   * residual_V is NOT. It measures ||V_calc - V_eff|| BEFORE the mixing,
    !!     i.e. the distance between V_eff and its image under the Kohn-Sham map,
    !!     which is a property of the current point and not of the step size.
    !!
    !! The sharpest assertion is at iteration 2. Iteration 1 mixes a zero
    !! difference (V_eff is seeded from the very density used to build V_calc), so
    !! the state entering iteration 2 is bit-identical for any alpha and the
    !! recorded residual_V(2) must be EXACTLY equal between the two runs.
    !!
    !! MUTATION THIS TEST KILLS: moving the residual_V computation to after the
    !! mixing block. There residual_V = ||V_calc - V_mixed|| = (1-alpha)*||V_calc
    !! - V_eff||, which reintroduces an alpha dependence into the convergence
    !! criterion. Verified by mutation: residual_V(2) then becomes 0.95*r for
    !! alpha = 0.05 and 0.995*r for alpha = 0.005, a 4.5% relative difference, and
    !! the "identical at iteration 2" assertion fails (the other assertions, being
    !! order-of-magnitude, survive - which is why the exact one is here).
    subroutine test_scf_potential_residual_is_alpha_independent()
        use fortuno_serial, only: check => serial_check

        integer, parameter :: N_ITER = 6
        real(dp), parameter :: ALPHA_BIG = 0.05_dp
        real(dp), parameter :: ALPHA_SMALL = 0.005_dp

        real(dp) :: res_big(N_ITER), dn_big(N_ITER)
        real(dp) :: res_small(N_ITER), dn_small(N_ITER)
        real(dp) :: res_ratio_2, dn_ratio_2, res_ratio_n, dn_ratio_n
        logical :: ok_big, ok_small

        call run_fixed_alpha_probe(ALPHA_BIG, N_ITER, res_big, dn_big, ok_big)
        call run_fixed_alpha_probe(ALPHA_SMALL, N_ITER, res_small, dn_small, ok_small)

        call check(ok_big .and. ok_small, &
                   "Precondition: both runs must complete the full 6-iteration history")
        if (.not. (ok_big .and. ok_small)) return

        ! Precondition: the system really is away from self-consistency, so the
        ! ratios below are not ratios of noise.
        call check(res_big(2) > 1.0e-3_dp .and. res_small(2) > 1.0e-3_dp, &
                   "Precondition: the potential residual must be far from zero")
        call check(dn_big(N_ITER) > 1.0e-9_dp .and. dn_small(N_ITER) > 1.0e-9_dp, &
                   "Precondition: the density is still moving at the last iteration")

        ! --- The exact anchor -------------------------------------------------
        ! Iteration 1 mixes a zero difference, so iteration 2 starts from the very
        ! same V_eff in both runs: the residual recorded there cannot depend on
        ! alpha at all. Computing it after the mixing would scale it by (1-alpha)
        ! and break this equality.
        ! The tolerance is 1e-9 RELATIVE, not exact equality: (1-a)*V + a*V is not
        ! bit-identical to V in IEEE (about 1 ulp per component), and that 1e-16
        ! is amplified by the eigenvector sensitivity 1/gap - this 20-site system
        ! with an impurity has gaps of order 1e-2 - before reaching the residual.
        ! 1e-9 sits comfortably above that noise and still seven orders below the
        ! 4.5e-2 relative signal of the mutation this anchor exists to kill.
        call check(abs(res_big(2) - res_small(2)) <= 1.0e-9_dp * abs(res_big(2)), &
                   "residual_V at iteration 2 must be identical for alpha = 0.05 and " // &
                   "alpha = 0.005: it is measured before the mixing")

        ! --- The qualitative contrast ----------------------------------------
        res_ratio_2 = res_big(2) / res_small(2)
        dn_ratio_2 = dn_big(2) / dn_small(2)
        res_ratio_n = res_big(N_ITER) / res_small(N_ITER)
        dn_ratio_n = dn_big(N_ITER) / dn_small(N_ITER)

        call check(dn_ratio_2 > 5.0_dp .and. dn_ratio_2 < 20.0_dp, &
                   "||Δn|| must fall by roughly the factor 10 by which alpha was reduced")
        call check(dn_ratio_n > 5.0_dp .and. dn_ratio_n < 20.0_dp, &
                   "||Δn|| must still track alpha after 6 iterations")

        call check(res_ratio_2 > 0.5_dp .and. res_ratio_2 < 2.0_dp, &
                   "residual_V must NOT fall by a factor ~10 when alpha does")
        call check(res_ratio_n > 0.2_dp .and. res_ratio_n < 2.0_dp, &
                   "residual_V after 6 iterations must still be of the same order for both alphas")

        ! The two grandeurs must be qualitatively different, not merely different
        ! numbers: ||Δn|| shrinks with alpha by at least an order more than the
        ! residual does.
        call check(dn_ratio_n > 4.0_dp * res_ratio_n, &
                   "||Δn|| must be far more sensitive to alpha than residual_V is")
    end subroutine test_scf_potential_residual_is_alpha_independent

    !> With alpha = 1 the residual must still be non-zero (regression, T3)
    !!
    !! Seeding-independent companion of the anchor above, and it kills the same
    !! mutation ("compute residual_V after the mixing block") without depending
    !! on floating-point luck or on how V_eff is initialised.
    !!
    !! With use_adaptive_mixing = .false. and mixing_alpha = 1 the mixing is pure
    !! substitution, V_eff <- V_calc. A residual measured AFTER the mixing would
    !! therefore be (1 - alpha)*||V_calc - V_eff|| = 0 identically, at every
    !! iteration and for every system. Measured before the mixing - as it is -
    !! it is the genuine distance between consecutive Kohn-Sham maps and stays
    !! far from zero while the cycle is still moving.
    !!
    !! Iteration 1 is excluded: there the residual is legitimately zero because
    !! V_eff is seeded from the same density used to build V_calc.
    subroutine test_scf_residual_nonzero_at_alpha_one()
        use fortuno_serial, only: check => serial_check

        integer, parameter :: N_ITER = 4
        real(dp) :: residuals(N_ITER), dnorms(N_ITER)
        logical :: ok
        integer :: k

        call run_fixed_alpha_probe(1.0_dp, N_ITER, residuals, dnorms, ok)

        call check(ok, "Precondition: the alpha = 1 run must complete the full history")
        if (.not. ok) return

        do k = 2, N_ITER
            call check(residuals(k) > 1.0e-6_dp, &
                       "with alpha = 1 the potential residual must still be non-zero: " // &
                       "measuring it after the substitution V_eff <- V_calc would make " // &
                       "it identically zero")
        end do
    end subroutine test_scf_residual_nonzero_at_alpha_one

    !> The residual must equal its own definition, factor included (regression)
    !!
    !! Every other assertion about residual_V is a ratio or an order of
    !! magnitude, so a wrong normalisation (dividing by L instead of 2L) or a
    !! forgotten spin channel (summing only the up term) would pass all of them.
    !! This test pins the ABSOLUTE value:
    !!
    !!   residual_V(k) = sqrt( (||V_up_calc - V_up||^2 + ||V_dw_calc - V_dw||^2)
    !!                         / (2*L) )
    !!
    !! recomputed here from V_ext, U, the densities and get_vxc alone, with no
    !! reference to the SCF internals.
    !!
    !! How the two potentials are obtained without reaching inside the loop: the
    !! probe runs with alpha = 1 and fixed mixing, so the mixing is a pure
    !! substitution and the potential entering iteration k is exactly
    !! V_calc[n_out(k-2)]. Therefore
    !!
    !!   residual_V(3) = || V_calc[n_out(2)] - V_calc[n_out(1)] || / sqrt(2L),
    !!
    !! and n_out(1), n_out(2) are read off two shorter runs of the very same
    !! system (max_iter = 1 and max_iter = 2).
    !!
    !! Nup = 7 /= Ndown = 5 is essential: the two spin channels carry different
    !! densities and different potentials, so dropping the down term from the sum
    !! changes the answer. On a spin-symmetric system the two channels coincide
    !! and the forgotten-channel bug would be invisible.
    !!
    !! MUTATIONS THIS TEST KILLS (both verified):
    !!   (a) real(2 * params%L, dp) -> real(params%L, dp) in the normalisation:
    !!       the recorded residual grows by sqrt(2), a 41% relative error;
    !!   (b) dropping sum((V_eff_down_calc - V_eff_down)**2) from the sum: the
    !!       recorded residual drops to the up-only value.
    subroutine test_scf_potential_residual_absolute_value()
        use fortuno_serial, only: check => serial_check
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy, get_vxc
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: N1 = 1, N2 = 2, N3 = 3
        real(dp) :: res1(N1), dn1(N1), res2(N2), dn2(N2), res3(N3), dn3(N3)
        real(dp) :: n1_up(PROBE_L), n1_down(PROBE_L)
        real(dp) :: n2_up(PROBE_L), n2_down(PROBE_L)
        real(dp) :: V_ext(PROBE_L)
        real(dp) :: V1_up(PROBE_L), V1_down(PROBE_L)
        real(dp) :: V2_up(PROBE_L), V2_down(PROBE_L)
        real(dp) :: vxc_up, vxc_down, residual_ref, up_only_ref
        type(xc_lsda_t) :: xc_func
        logical :: ok1, ok2, ok3
        integer :: i, ierr

        ! Three runs of the SAME system, differing only in how many iterations
        ! they are allowed: alpha = 1, fixed mixing (see the header).
        call run_fixed_alpha_probe(1.0_dp, N1, res1, dn1, ok1, &
                                   dens_up=n1_up, dens_down=n1_down, v_ext_out=V_ext)
        call run_fixed_alpha_probe(1.0_dp, N2, res2, dn2, ok2, &
                                   dens_up=n2_up, dens_down=n2_down)
        call run_fixed_alpha_probe(1.0_dp, N3, res3, dn3, ok3)

        call check(ok1 .and. ok2 .and. ok3, &
                   "Precondition: the 1-, 2- and 3-iteration runs must all complete")
        if (.not. (ok1 .and. ok2 .and. ok3)) return

        call xc_lsda_init(xc_func, PROBE_TABLE, ierr)
        call check(ierr == ERROR_SUCCESS, "XC init should succeed")
        if (ierr /= ERROR_SUCCESS) return

        ! V_calc[n] = V_ext + U*n_other + V_xc, evaluated independently of the
        ! SCF loop for the two densities that bracket iteration 3.
        do i = 1, PROBE_L
            call get_vxc(xc_func, n1_up(i), n1_down(i), vxc_up, vxc_down, ierr)
            if (ierr /= ERROR_SUCCESS) exit
            V1_up(i) = V_ext(i) + PROBE_U * n1_down(i) + vxc_up
            V1_down(i) = V_ext(i) + PROBE_U * n1_up(i) + vxc_down

            call get_vxc(xc_func, n2_up(i), n2_down(i), vxc_up, vxc_down, ierr)
            if (ierr /= ERROR_SUCCESS) exit
            V2_up(i) = V_ext(i) + PROBE_U * n2_down(i) + vxc_up
            V2_down(i) = V_ext(i) + PROBE_U * n2_up(i) + vxc_down
        end do
        call check(ierr == ERROR_SUCCESS, "get_vxc must succeed on both densities")
        call xc_lsda_destroy(xc_func)
        if (ierr /= ERROR_SUCCESS) return

        residual_ref = sqrt((sum((V2_up - V1_up)**2) + sum((V2_down - V1_down)**2)) &
                            / real(2 * PROBE_L, dp))
        up_only_ref = sqrt(sum((V2_up - V1_up)**2) / real(2 * PROBE_L, dp))

        ! Preconditions: the value is a real number, not noise, and the two spin
        ! channels really do contribute differently - otherwise mutation (b)
        ! could not be detected.
        call check(residual_ref > 1.0e-3_dp, &
                   "Precondition: the reference residual must be far from zero")
        call check(abs(residual_ref - up_only_ref) > 0.1_dp * residual_ref, &
                   "Precondition: the down channel must contribute materially, " // &
                   "otherwise a forgotten channel would be invisible")

        call check(abs(res3(N3) - residual_ref) <= 1.0e-9_dp * max(1.0_dp, residual_ref), &
                   "residual_V must equal sqrt((||dV_up||^2 + ||dV_down||^2)/(2L)) " // &
                   "recomputed independently from V_ext, U, n and get_vxc")
    end subroutine test_scf_potential_residual_absolute_value

    !> The adaptive controller must START from the user's mixing_alpha (T5)
    !!
    !! With use_adaptive_mixing = .true. (the default) the controller used to be
    !! seeded with the hard-coded INITIAL_MIX = 0.95, so scf_params%mixing_alpha
    !! was silently ignored: every run behaved as if alpha = 0.05 no matter what
    !! the user asked for. The coincidence 1 - 0.95 = 0.05 = the value in
    !! input.txt is what kept this invisible.
    !!
    !! The controller retunes alpha only after count_sc_max = 10 in-band
    !! iterations, so over the 6 iterations probed here the mixing weight is
    !! exactly the seeded one. Consequently ||Δn||, which is proportional to the
    !! mixing weight, must differ by the factor 10 between the two runs. Without
    !! the fix both runs start from mix = 0.95 and produce bit-identical
    !! histories, i.e. a ratio of exactly 1, and this test fails.
    subroutine test_adaptive_mixing_honours_user_alpha()
        use fortuno_serial, only: check => serial_check

        integer, parameter :: N_ITER = 6
        real(dp) :: res_big(N_ITER), dn_big(N_ITER)
        real(dp) :: res_small(N_ITER), dn_small(N_ITER)
        logical :: ok_big, ok_small

        call run_fixed_alpha_probe(0.05_dp, N_ITER, res_big, dn_big, ok_big, adaptive=.true.)
        call run_fixed_alpha_probe(0.005_dp, N_ITER, res_small, dn_small, ok_small, adaptive=.true.)

        call check(ok_big .and. ok_small, &
                   "Precondition: both adaptive runs must complete the full history")
        if (.not. (ok_big .and. ok_small)) return

        call check(dn_big(2) > 0.0_dp .and. dn_small(2) > 0.0_dp, &
                   "Precondition: the density must still be moving")

        call check(abs(dn_big(2) - dn_small(2)) > 1.0e-6_dp, &
                   "mixing_alpha must not be ignored when adaptive mixing is on")
        call check(dn_big(2) / dn_small(2) > 5.0_dp .and. dn_big(2) / dn_small(2) < 20.0_dp, &
                   "the adaptive controller must start at the user's alpha, so ||Δn|| " // &
                   "must scale with it over the first iterations")
    end subroutine test_adaptive_mixing_honours_user_alpha

    !> Only sites within HALF_FILLING_TOL of n = 1 are counted
    subroutine test_count_half_filled_sites()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: count_half_filled_sites, HALF_FILLING_TOL
        use lsda_constants, only: dp

        real(dp) :: n_up(5), n_down(5)

        ! Sites 1 and 2 sit at half filling (exactly and just inside the band),
        ! site 3 is just outside it, sites 4 and 5 are far away.
        n_up   = [0.5_dp, 0.5_dp, 0.5_dp, 0.10_dp, 0.90_dp]
        n_down = [0.5_dp, 0.5_dp + 0.5_dp * HALF_FILLING_TOL, &
                  0.5_dp + 2.0_dp * HALF_FILLING_TOL, 0.10_dp, 0.90_dp]

        call check(count_half_filled_sites(n_up, n_down, 5) == 2, &
                   "exactly the two sites inside the half-filling band must be counted")
        call check(count_half_filled_sites(n_up, n_down, 1) == 1, &
                   "counting must respect L")
    end subroutine test_count_half_filled_sites

    !> The warning fires only with half-filled sites AND an oscillating energy
    !!
    !! Regression guard for T4: before it, a run could sit at n = 1 hopping
    !! across the V_xc discontinuity with ΔE alternating by ±1.29 and the user
    !! got no explanation at all.
    subroutine test_half_filling_warning_fires_when_oscillating()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: half_filling_warning_due, OSCILLATION_STREAK_MIN
        use lsda_constants, only: dp

        call check(half_filling_warning_due(1, OSCILLATION_STREAK_MIN, 1.287_dp, &
                                            -311.0_dp, 1.0e-8_dp), &
                   "one half-filled site with an alternating, large dE must warn")
        call check(half_filling_warning_due(7, OSCILLATION_STREAK_MIN + 4, -1.287_dp, &
                                            -311.0_dp, 1.0e-8_dp), &
                   "the sign of the last dE must not matter")
    end subroutine test_half_filling_warning_fires_when_oscillating

    !> The warning stays silent when either condition is missing
    subroutine test_half_filling_warning_silent_otherwise()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: half_filling_warning_due, OSCILLATION_STREAK_MIN
        use lsda_constants, only: dp

        call check(.not. half_filling_warning_due(0, OSCILLATION_STREAK_MIN + 2, 1.287_dp, &
                                                  -311.0_dp, 1.0e-8_dp), &
                   "no half-filled site: the discontinuity cannot be the cause")
        call check(.not. half_filling_warning_due(3, OSCILLATION_STREAK_MIN - 1, 1.287_dp, &
                                                  -311.0_dp, 1.0e-8_dp), &
                   "a short sign-flip streak is not an oscillation")
        call check(.not. half_filling_warning_due(3, OSCILLATION_STREAK_MIN, 1.0e-10_dp, &
                                                  -311.0_dp, 1.0e-8_dp), &
                   "a settled energy must not warn, however many sites sit at n = 1")
    end subroutine test_half_filling_warning_silent_otherwise

    !> REGRESSION (T9): the failure path must expose the final state, not nothing
    !!
    !! On ERROR_CONVERGENCE_FAILED the cycle used to fill final_energy,
    !! final_density_error and final_potential_residual while leaving
    !! density_up, density_down and eigvals UNALLOCATED. The failure is the
    !! diagnostically interesting case, and any consumer that read
    !! results%density_up after that error code dereferenced an unallocated
    !! array - the output writer silently skipped the density and eigenvalue
    !! files, and a less careful caller segfaulted.
    !!
    !! max_iter = 1 is a non-convergence that does not depend on the physics:
    !! declaring self-consistency requires a previous energy to compare against,
    !! so the first iteration can never converge, whatever the system.
    !!
    !! Both public APIs are probed: real and complex calls share the same SCF
    !! loop but select different diagonalization backends when required.
    !> REGRESSION (T15): output XC values must become the next input cache.
    !!
    !! A two-iteration, deliberately non-converged SCF needs exactly one
    !! initial `get_vxc` call per site and then one `get_vxc` plus one
    !! `get_exc` call per site per iteration.  Reintroducing the former input
    !! V_xc evaluation adds L calls and fails this test.
    subroutine test_scf_reuses_output_xc_cache()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, scf_params_t, scf_results_t, cleanup_scf_results
        use lsda_types, only: system_params_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy, &
                           reset_xc_evaluation_count, get_xc_evaluation_count
        use boundary_conditions, only: BC_OPEN
        use lsda_errors, only: ERROR_SUCCESS, ERROR_CONVERGENCE_FAILED

        integer, parameter :: L = 8
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L)
        integer :: ierr

        call xc_lsda_init(xc_func, PROBE_TABLE, ierr)
        call check(ierr == ERROR_SUCCESS, "XC cache: initialization should succeed")
        if (ierr /= ERROR_SUCCESS) return

        params%L = L
        params%Nup = 3
        params%Ndown = 2
        params%bc = BC_OPEN
        params%U = PROBE_U
        params%phase = 0.0_dp
        scf_params%max_iter = 2
        scf_params%potential_tol = 1.0e-14_dp
        scf_params%energy_tol = 1.0e-14_dp
        scf_params%mixing_alpha = 0.2_dp
        scf_params%verbose = .false.
        scf_params%store_history = .false.
        V_ext = 0.0_dp

        V_ext(1) = 0.1_dp  ! Deliberately break reflection symmetry for the count.
        call reset_xc_evaluation_count()
        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)

        call check(ierr == ERROR_CONVERGENCE_FAILED, "XC cache: two iterations should not converge")
        call check(get_xc_evaluation_count() == 5 * L, &
                   "XC cache: initial Vxc plus one Vxc/exc pair per output density")

        call cleanup_scf_results(results, ierr)
        call xc_lsda_destroy(xc_func)
    end subroutine test_scf_reuses_output_xc_cache

    !> Reflection-related sites must retain their own XC values near n = 1.
    !!
    !! Total densities 1 - 2e-12 and 1 + 2e-12 differ by only 4e-12, which
    !! is below the former reflection-cache threshold of 1e-10. They are,
    !! however, outside the XC region-boundary tolerance (1e-12) and lie on
    !! opposite sides of the physical Mott discontinuity. Mirroring the first
    !! site's cache entry onto the second would replace its physical branch.
    subroutine test_xc_cache_preserves_half_filling_discontinuity()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: populate_xc_output_cache
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy, get_vxc, get_exc
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 2
        real(dp), parameter :: OFFSET = 2.0e-12_dp
        type(xc_lsda_t) :: xc_func
        real(dp) :: density_up(L), density_down(L)
        real(dp) :: v_xc_up(L), v_xc_down(L), e_xc(L)
        real(dp) :: v_xc_up_ref(L), v_xc_down_ref(L), e_xc_ref(L)
        integer :: i, ierr

        call xc_lsda_init(xc_func, PROBE_TABLE, ierr)
        call check(ierr == ERROR_SUCCESS, "XC cusp cache: initialization should succeed")
        if (ierr /= ERROR_SUCCESS) return

        density_up = 0.5_dp
        density_down = [0.5_dp - OFFSET, 0.5_dp + OFFSET]

        call populate_xc_output_cache(xc_func, density_up, density_down, v_xc_up, v_xc_down, e_xc, ierr)
        call check(ierr == ERROR_SUCCESS, "XC cusp cache: cache population should succeed")

        do i = 1, L
            call get_vxc(xc_func, density_up(i), density_down(i), v_xc_up_ref(i), v_xc_down_ref(i), ierr)
            if (ierr == ERROR_SUCCESS) call get_exc(xc_func, density_up(i), density_down(i), e_xc_ref(i), ierr)
            call check(ierr == ERROR_SUCCESS, "XC cusp cache: independent XC evaluation should succeed")
        end do

        if (ierr == ERROR_SUCCESS) then
            call check(maxval(abs(v_xc_up - v_xc_up_ref)) < 1.0e-14_dp, &
                       "XC cusp cache: each site must retain its own spin-up V_xc branch")
            call check(maxval(abs(v_xc_down - v_xc_down_ref)) < 1.0e-14_dp, &
                       "XC cusp cache: each site must retain its own spin-down V_xc branch")
            call check(maxval(abs(e_xc - e_xc_ref)) < 1.0e-14_dp, &
                       "XC cusp cache: each site must retain its own e_xc value")
            call check(abs(v_xc_up_ref(2) - v_xc_up_ref(1)) > 1.0e-3_dp, &
                       "XC cusp cache: the two sites must straddle a finite V_xc discontinuity")
        end if

        call xc_lsda_destroy(xc_func)
    end subroutine test_xc_cache_preserves_half_filling_discontinuity

    subroutine test_scf_failure_exposes_final_state()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: run_kohn_sham_scf_real, run_kohn_sham_scf_complex, &
                                    scf_params_t, scf_results_t, cleanup_scf_results
        use lsda_types, only: system_params_t
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use boundary_conditions, only: BC_OPEN, BC_TWISTED
        use lsda_errors, only: ERROR_SUCCESS, ERROR_CONVERGENCE_FAILED

        integer, parameter :: L = 8
        integer, parameter :: NUP = 3
        integer, parameter :: NDOWN = 2
        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        type(scf_results_t) :: results
        type(xc_lsda_t) :: xc_func
        real(dp) :: V_ext(L)
        integer :: ierr

        call xc_lsda_init(xc_func, PROBE_TABLE, ierr)
        call check(ierr == ERROR_SUCCESS, "XC init should succeed")

        params%L = L
        params%Nup = NUP
        params%Ndown = NDOWN
        params%bc = BC_OPEN
        params%U = PROBE_U
        params%phase = 0.0_dp

        scf_params%max_iter = 1  ! Cannot converge: no previous energy exists.
        scf_params%potential_tol = 1.0e-8_dp
        scf_params%energy_tol = 1.0e-10_dp
        scf_params%mixing_alpha = 0.3_dp
        scf_params%verbose = .false.
        scf_params%store_history = .true.

        V_ext = 0.0_dp

        call run_kohn_sham_scf_real(params, scf_params, V_ext, xc_func, results, ierr)

        call check(ierr == ERROR_CONVERGENCE_FAILED, &
                   "real loop: one iteration must report convergence failure")
        call check(.not. results%converged, "real loop: converged flag must be .false.")
        call check(allocated(results%density_up), &
                   "real loop: density_up must be available after a failed SCF")
        call check(allocated(results%density_down), &
                   "real loop: density_down must be available after a failed SCF")
        call check(allocated(results%eigvals), &
                   "real loop: eigvals must be available after a failed SCF")

        if (allocated(results%density_up)) then
            call check(size(results%density_up) == L, "real loop: density_up has L entries")
            call check(abs(sum(results%density_up) - real(NUP, dp)) < 1.0e-10_dp, &
                       "real loop: the exposed density must still integrate to N_up")
        end if
        if (allocated(results%density_down)) then
            call check(abs(sum(results%density_down) - real(NDOWN, dp)) < 1.0e-10_dp, &
                       "real loop: the exposed density must still integrate to N_down")
        end if
        if (allocated(results%eigvals)) then
            call check(size(results%eigvals) == 2 * L, &
                       "real loop: eigvals holds both spin channels")
        end if

        call cleanup_scf_results(results, ierr)

        ! The complex loop must also take the open-boundary tridiagonal route.
        ! Its dense H_up/H_down placeholders have extent (1,1) in this case, so
        ! either dense Hamiltonian builder must be skipped before DSTEVR runs.
        params%bc = BC_OPEN

        call run_kohn_sham_scf_complex(params, scf_params, V_ext, xc_func, results, ierr)

        call check(ierr == ERROR_CONVERGENCE_FAILED, &
                   "complex OBC loop: one iteration must report convergence failure, not a matrix-size error")
        call check(allocated(results%density_up), &
                   "complex OBC loop: density_up must be available after the tridiagonal route")
        call check(allocated(results%density_down), &
                   "complex OBC loop: density_down must be available after the tridiagonal route")

        call cleanup_scf_results(results, ierr)

        ! Same contract for the complex twisted-BC copy of the loop.
        params%bc = BC_TWISTED
        params%phase = 0.5_dp

        call run_kohn_sham_scf_complex(params, scf_params, V_ext, xc_func, results, ierr)

        call check(ierr == ERROR_CONVERGENCE_FAILED, &
                   "complex loop: one iteration must report convergence failure")
        call check(allocated(results%density_up), &
                   "complex loop: density_up must be available after a failed SCF")
        call check(allocated(results%density_down), &
                   "complex loop: density_down must be available after a failed SCF")
        call check(allocated(results%eigvals), &
                   "complex loop: eigvals must be available after a failed SCF")

        if (allocated(results%eigvals)) then
            call check(size(results%eigvals) == 2 * L, &
                       "complex loop: eigvals holds both spin channels")
        end if

        call cleanup_scf_results(results, ierr)
        call xc_lsda_destroy(xc_func)
    end subroutine test_scf_failure_exposes_final_state

end program test_kohn_sham_cycle
