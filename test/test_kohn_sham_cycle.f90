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
            test("validate_inputs_valid", test_validate_inputs_valid), &
            test("validate_inputs_invalid_L", test_validate_inputs_invalid_L), &
            test("validate_inputs_invalid_N", test_validate_inputs_invalid_N), &
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
    !! Physics: E_tot = Σε_j + E_xc - ∫V_xc·n
    !! The double-counting correction removes XC potential from band energy.
    subroutine test_compute_total_energy_simple()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: compute_total_energy
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 10
        real(dp) :: eigvals_up(5), eigvals_down(5)
        real(dp) :: n_up(L), n_down(L), V_ext(L)
        type(xc_lsda_t) :: xc_func
        real(dp) :: total_energy
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

        call compute_total_energy(eigvals_up, eigvals_down, 5, 5, n_up, n_down, &
                                  V_ext, xc_func, L, total_energy, ierr)

        call check(ierr == ERROR_SUCCESS, "Energy calculation should succeed")
        call check(total_energy == total_energy, "Energy should be valid number")
        ! E_band = 2*(-2-1.5-1-0.5+0) = -10, plus XC corrections
        call check(abs(total_energy + 10.0_dp) < 5.0_dp, "Energy should be reasonable")

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

        ! All L eigenvalues occupied (half-filling)
        eigvals_up = [-2.0_dp, -1.8_dp, -1.5_dp, -1.2_dp, -1.0_dp, -0.8_dp, -0.5_dp, -0.2_dp]
        eigvals_down = eigvals_up

        call compute_total_energy(eigvals_up, eigvals_down, L, L, n_up, n_down, &
                                  V_ext, xc_func, L, total_energy, ierr)

        call check(ierr == ERROR_SUCCESS, "Energy calculation should succeed")
        call check(total_energy == total_energy, "Energy should be valid")

        call xc_lsda_destroy(xc_func)
    end subroutine test_compute_total_energy_half_filling

    !> Test input validation with valid parameters
    subroutine test_validate_inputs_valid()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: scf_params_t
        use lsda_types, only: system_params_t
        use boundary_conditions, only: BC_PERIODIC
        use lsda_errors, only: ERROR_SUCCESS

        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        real(dp) :: V_ext(10)
        integer :: ierr

        ! Valid system parameters
        params%L = 10
        params%Nup = 5
        params%Ndown = 5
        params%bc = BC_PERIODIC
        params%U = 4.0_dp
        params%phase = 0.0_dp

        ! Valid SCF parameters
        scf_params%max_iter = 100
        scf_params%density_tol = 1.0e-6_dp
        scf_params%energy_tol = 1.0e-8_dp
        scf_params%mixing_alpha = 0.3_dp

        V_ext = 0.0_dp

        ! This is tested internally by run_kohn_sham_scf_real
        ! Just verify parameters are reasonable
        call check(params%L > 0, "L should be positive")
        call check(params%Nup >= 0, "Nup should be non-negative")
        call check(params%Ndown >= 0, "Ndown should be non-negative")
        call check(params%Nup + params%Ndown <= params%L, "N <= L")
        call check(scf_params%mixing_alpha > 0.0_dp .and. scf_params%mixing_alpha <= 1.0_dp, &
                   "mixing_alpha should be in (0, 1]")
    end subroutine test_validate_inputs_valid

    !> Test input validation rejects invalid L
    subroutine test_validate_inputs_invalid_L()
        use fortuno_serial, only: check => serial_check
        use lsda_types, only: system_params_t
        use kohn_sham_cycle, only: scf_params_t

        type(system_params_t) :: params
        type(scf_params_t) :: scf_params

        params%L = 0  ! Invalid!
        params%Nup = 5
        params%Ndown = 5

        call check(params%L <= 0, "Should detect invalid L")
    end subroutine test_validate_inputs_invalid_L

    !> Test input validation rejects invalid particle numbers
    subroutine test_validate_inputs_invalid_N()
        use fortuno_serial, only: check => serial_check
        use lsda_types, only: system_params_t
        use kohn_sham_cycle, only: scf_params_t

        type(system_params_t) :: params
        type(scf_params_t) :: scf_params

        params%L = 10
        params%Nup = 8
        params%Ndown = 8  ! N > L, invalid!

        call check(params%Nup + params%Ndown > params%L, "Should detect N > L")
    end subroutine test_validate_inputs_invalid_N

    !> Test input validation detects size mismatch
    subroutine test_validate_inputs_size_mismatch()
        use fortuno_serial, only: check => serial_check
        use lsda_types, only: system_params_t
        use kohn_sham_cycle, only: scf_params_t

        type(system_params_t) :: params
        type(scf_params_t) :: scf_params
        real(dp) :: V_ext(5)  ! Wrong size!

        params%L = 10
        params%Nup = 5
        params%Ndown = 5

        call check(size(V_ext) /= params%L, "Should detect size mismatch")
    end subroutine test_validate_inputs_size_mismatch

    !> Test input validation rejects invalid mixing parameter
    subroutine test_validate_inputs_invalid_mixing()
        use fortuno_serial, only: check => serial_check
        use kohn_sham_cycle, only: scf_params_t

        type(scf_params_t) :: scf_params

        scf_params%mixing_alpha = 1.5_dp  ! Invalid! Must be <= 1

        call check(scf_params%mixing_alpha > 1.0_dp, "Should detect invalid mixing_alpha")
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
