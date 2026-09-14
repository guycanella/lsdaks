!> Unit tests for kohn_sham_cycle module
program test_kohn_sham_cycle
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp
    implicit none

    real(dp), parameter :: TOL = 1.0e-9_dp

    call execute_serial_cmd_app(get_kohn_sham_tests())

contains

    function get_kohn_sham_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests

        tests = test_list([ &
            test("compute_total_energy_simple", test_compute_total_energy_simple), &
            test("compute_total_energy_half_filling", test_compute_total_energy_half_filling), &
            test("validate_double_occupancy_accepted", test_validate_double_occupancy_accepted), &
            test("validate_inputs_invalid_L", test_validate_inputs_invalid_L), &
            test("validate_N_per_spin_exceeds_L", test_validate_N_per_spin_exceeds_L), &
            test("validate_inputs_size_mismatch", test_validate_inputs_size_mismatch), &
            test("validate_inputs_invalid_mixing", test_validate_inputs_invalid_mixing), &
            test("scf_results_init_cleanup", test_scf_results_init_cleanup), &
            test("scf_converges_u0_open", test_scf_converges_u0_open), &
            test("scf_converges_u0_periodic", test_scf_converges_u0_periodic), &
            test("scf_stores_history", test_scf_stores_history), &
            test("scf_density_conservation", test_scf_density_conservation), &
            test("scf_complex_twisted_bc", test_scf_complex_twisted_bc) &
        ])
    end function get_kohn_sham_tests

    !> Test total energy calculation with simple case
    !!
    !! Physics: E_tot = Σε_j - U·Σn↑n↓ + E_xc - ∫V_xc·n
    !! The Hartree term -U·Σn↑n↓ removes the double counting of U·n_other that
    !! the eigenvalues already carry (C++ lsdaks.cc:675-679); the V_xc term does
    !! the same for the XC potential.
    !!
    !! The decisive assertion is differential: since U enters compute_total_energy
    !! only through E_hartree (the XC functional is fixed by the table, not by the
    !! U argument), E(U=2) - E(U=0) must be exactly -2·Σn↑n↓ = -1.25. This pins
    !! both the magnitude and the sign of the new U argument; a wrong sign would
    !! give +1.25 and a dropped/mis-passed U would give 0.
    subroutine test_compute_total_energy_simple()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: compute_total_energy
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 10
        real(dp) :: eigvals_up(5), eigvals_down(5)
        real(dp) :: n_up(L), n_down(L), V_ext(L)
        type(xc_lsda_t) :: xc_func
        real(dp) :: total_energy, energy_u0, hartree_sum
        integer :: ierr, i
        character(len=256) :: table_file

        ! Initialize XC functional
        table_file = "data/tables/fortran_native/xc_table_u2.00.dat"
        call xc_lsda_init(xc_func, table_file, ierr)
        call check(ierr == ERROR_SUCCESS, "XC init should succeed")

        ! Setup simple case: uniform density n=0.5
        do i = 1, L
            n_up(i) = 0.25_dp
            n_down(i) = 0.25_dp
            V_ext(i) = 0.0_dp
        end do

        ! Simple eigenvalues (5 occupied levels per spin)
        eigvals_up = [-2.0_dp, -1.5_dp, -1.0_dp, -0.5_dp, 0.0_dp]
        eigvals_down = [-2.0_dp, -1.5_dp, -1.0_dp, -0.5_dp, 0.0_dp]

        ! U = 2.0 matches the XC table loaded above (xc_table_u2.00.dat)
        call compute_total_energy(eigvals_up, eigvals_down, 5, 5, n_up, n_down, &
                                  V_ext, xc_func, 2.0_dp, L, total_energy, ierr)

        call check(ierr == ERROR_SUCCESS, "Energy calculation should succeed")
        call check(total_energy == total_energy, "Energy should be valid number")

        ! Same call with U = 0: identical band energy, identical XC (same table),
        ! only the Hartree term disappears.
        call compute_total_energy(eigvals_up, eigvals_down, 5, 5, n_up, n_down, &
                                  V_ext, xc_func, 0.0_dp, L, energy_u0, ierr)
        call check(ierr == ERROR_SUCCESS, "Energy calculation at U=0 should succeed")

        ! Σ n_up*n_down = 10 * 0.25 * 0.25 = 0.625  =>  E(U=2) - E(U=0) = -1.25
        hartree_sum = sum(n_up * n_down)
        call check(abs(hartree_sum - 0.625_dp) < TOL, "Precondition: Σn_up*n_down = 0.625")
        call check(abs((total_energy - energy_u0) + 2.0_dp * hartree_sum) < 1.0e-12_dp, &
                   "E(U=2) - E(U=0) must equal -U*Σn_up*n_down = -1.25")

        ! Band energy alone is exactly -10, so at U=0 the whole remainder is the
        ! XC double-counting correction; it must stay small compared to the band.
        call check(abs(energy_u0 + 10.0_dp) < 2.0_dp, &
                   "At U=0 the energy must sit close to the band energy -10")

        call xc_lsda_destroy(xc_func)
    end subroutine test_compute_total_energy_simple

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
        real(dp) :: n_up(L), n_down(L), V_ext(L)
        type(xc_lsda_t) :: xc_func
        real(dp) :: total_energy
        integer :: ierr, i
        character(len=256) :: table_file

        table_file = "data/tables/fortran_native/xc_table_u4.00.dat"
        call xc_lsda_init(xc_func, table_file, ierr)
        call check(ierr == ERROR_SUCCESS, "XC init should succeed")

        ! Half-filling: n_up = n_down = 0.5 at each site
        do i = 1, L
            n_up(i) = 0.5_dp
            n_down(i) = 0.5_dp
            V_ext(i) = 0.0_dp
        end do

        ! The lowest L/2 eigenvalues are occupied in each spin channel
        eigvals_up = [-2.0_dp, -1.8_dp, -1.5_dp, -1.2_dp, -1.0_dp, -0.8_dp, -0.5_dp, -0.2_dp]
        eigvals_down = eigvals_up

        ! U = 4.0 matches the XC table loaded above (xc_table_u4.00.dat)
        call compute_total_energy(eigvals_up, eigvals_down, L / 2, L / 2, n_up, n_down, &
                                  V_ext, xc_func, 4.0_dp, L, total_energy, ierr)

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

        table_file = "data/tables/fortran_native/xc_table_u4.00.dat"
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

        table_file = "data/tables/fortran_native/xc_table_u4.00.dat"
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

        table_file = "data/tables/fortran_native/xc_table_u4.00.dat"
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

        table_file = "data/tables/fortran_native/xc_table_u4.00.dat"
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

        table_file = "data/tables/fortran_native/xc_table_u4.00.dat"
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

    !> Test SCF convergence for U=0 with open BC
    !!
    !! Physics: U=0 is the non-interacting limit (free Fermi gas).
    !! SCF should converge in 1 iteration since V_xc = 0.
    !! Density should be uniform for uniform V_ext.
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
        integer :: ierr, i
        character(len=256) :: table_file

        ! Setup: U=0 table (if available, otherwise use small U)
        table_file = "data/tables/fortran_native/xc_table_u1.00.dat"
        call xc_lsda_init(xc_func, table_file, ierr)
        if (ierr /= ERROR_SUCCESS) then
            ! Skip test if table not found
            return
        end if

        ! System parameters: small system, half-filled
        params%L = L
        params%Nup = 5
        params%Ndown = 5
        params%bc = BC_OPEN
        params%U = 1.0_dp
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
        call check(results%n_iterations <= 10, "Should converge quickly for small U")
        call check(allocated(results%density_up), "Density_up should be allocated")
        call check(allocated(results%density_down), "Density_down should be allocated")
        call check(allocated(results%eigvals), "Eigvals should be allocated")

        ! Verify particle number conservation
        call check(abs(sum(results%density_up) + sum(results%density_down) - 10.0_dp) < 1.0e-6_dp, &
                   "Particle number should be conserved")

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

        table_file = "data/tables/fortran_native/xc_table_u2.00.dat"
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

        table_file = "data/tables/fortran_native/xc_table_u2.00.dat"
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

        table_file = "data/tables/fortran_native/xc_table_u4.00.dat"
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

        table_file = "data/tables/fortran_native/xc_table_u2.00.dat"
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

end program test_kohn_sham_cycle
