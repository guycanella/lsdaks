!> Unit tests for convergence_monitor module
program test_convergence_monitor
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp
    implicit none

    real(dp), parameter :: TOL = 1.0e-9_dp

    call execute_serial_cmd_app(get_convergence_tests())

contains

    function get_convergence_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests

        tests = test_list([ &
            test("density_difference_simple", test_density_difference_simple), &
            test("density_difference_zero", test_density_difference_zero), &
            test("density_difference_size_mismatch", test_density_difference_size_mismatch), &
            test("density_norm_L1", test_density_norm_L1), &
            test("density_norm_L2", test_density_norm_L2), &
            test("density_norm_Linf", test_density_norm_Linf), &
            test("density_norm_invalid_type", test_density_norm_invalid_type), &
            test("convergence_check_converged", test_convergence_check_converged), &
            test("convergence_check_not_converged", test_convergence_check_not_converged), &
            test("convergence_check_custom_tolerance", test_convergence_check_custom_tolerance), &
            test("history_init_cleanup", test_history_init_cleanup), &
            test("history_update", test_history_update), &
            test("history_bounds_checking", test_history_bounds_checking) &
        ])
    end function get_convergence_tests

    !> Test simple density difference calculation
    !!
    !! Physics: During SCF iteration, we need to compute Δn = n_out - n_in
    !! to monitor convergence. This is the most basic operation in SCF cycles.
    subroutine test_density_difference_simple()
        use fortuno_serial, only: check => serial_check
        use convergence_monitor, only: compute_density_difference
        integer, parameter :: L = 5
        real(dp) :: n_new(L), n_old(L), diff(L)
        integer :: ierr

        ! Setup: new density slightly higher than old
        n_old = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
        n_new = [1.1_dp, 0.9_dp, 1.05_dp, 0.95_dp, 1.0_dp]

        call compute_density_difference(n_new, n_old, L, diff, ierr)

        call check(ierr == 0, "Difference: computation should succeed")
        call check(abs(diff(1) - 0.1_dp) < TOL, "Difference: site 1 should be 0.1")
        call check(abs(diff(2) - (-0.1_dp)) < TOL, "Difference: site 2 should be -0.1")
        call check(abs(diff(3) - 0.05_dp) < TOL, "Difference: site 3 should be 0.05")
        call check(abs(diff(5)) < TOL, "Difference: site 5 should be 0 (no change)")
    end subroutine test_density_difference_simple

    !> Test density difference when densities are identical
    !!
    !! Physics: When SCF has converged, n_in = n_out, so Δn = 0 everywhere.
    !! This is the convergence criterion.
    subroutine test_density_difference_zero()
        use fortuno_serial, only: check => serial_check
        use convergence_monitor, only: compute_density_difference
        integer, parameter :: L = 4
        real(dp) :: n_new(L), n_old(L), diff(L)
        integer :: ierr, i

        ! Identical densities (converged state)
        n_old = [0.8_dp, 1.2_dp, 0.5_dp, 1.0_dp]
        n_new = n_old

        call compute_density_difference(n_new, n_old, L, diff, ierr)

        call check(ierr == 0, "Zero diff: computation should succeed")
        do i = 1, L
            call check(abs(diff(i)) < TOL, "Zero diff: all differences should be zero")
        end do
    end subroutine test_density_difference_zero

    !> Test size mismatch detection
    !!
    !! Physics: Arrays must match system size L. Mismatches indicate
    !! programming errors or corrupted data structures.
    subroutine test_density_difference_size_mismatch()
        use fortuno_serial, only: check => serial_check
        use convergence_monitor, only: compute_density_difference
        use lsda_errors, only: ERROR_SIZE_MISMATCH
        integer, parameter :: L = 5
        real(dp) :: n_new(L), n_old(L), diff(L+1)  ! Wrong size!
        integer :: ierr

        n_new = 1.0_dp
        n_old = 1.0_dp

        call compute_density_difference(n_new, n_old, L, diff, ierr)

        call check(ierr == ERROR_SIZE_MISMATCH, &
                   "Size mismatch: should return ERROR_SIZE_MISMATCH")
    end subroutine test_density_difference_size_mismatch

    !> Test L1 norm calculation
    !!
    !! Physics: L1 norm = Σ|Δn(i)| measures total absolute density change.
    !! Useful for detecting charge sloshing in SCF.
    subroutine test_density_norm_L1()
        use fortuno_serial, only: check => serial_check
        use convergence_monitor, only: compute_density_norm, L1
        integer, parameter :: L = 4
        real(dp) :: delta_n(L), norm_value
        real(dp) :: expected_L1
        integer :: ierr

        ! Example: Δn = [0.1, -0.2, 0.3, -0.1]
        delta_n = [0.1_dp, -0.2_dp, 0.3_dp, -0.1_dp]
        expected_L1 = 0.1_dp + 0.2_dp + 0.3_dp + 0.1_dp  ! = 0.7

        call compute_density_norm(delta_n, L, L1, norm_value, ierr)

        call check(ierr == 0, "L1 norm: computation should succeed")
        call check(abs(norm_value - expected_L1) < TOL, "L1 norm: should equal 0.7")
    end subroutine test_density_norm_L1

    !> Test L2 norm calculation
    !!
    !! Physics: L2 norm = √(Σ|Δn(i)|²) is the Euclidean distance between densities.
    !! Standard convergence criterion in DFT: ||Δn||₂ < tolerance.
    subroutine test_density_norm_L2()
        use fortuno_serial, only: check => serial_check
        use convergence_monitor, only: compute_density_norm, L2
        integer, parameter :: L = 3
        real(dp) :: delta_n(L), norm_value
        real(dp) :: expected_L2
        integer :: ierr

        ! Example: Δn = [0.3, 0.4, 0.0]
        delta_n = [0.3_dp, 0.4_dp, 0.0_dp]
        expected_L2 = sqrt(0.3_dp**2 + 0.4_dp**2)  ! = sqrt(0.25) = 0.5

        call compute_density_norm(delta_n, L, L2, norm_value, ierr)

        call check(ierr == 0, "L2 norm: computation should succeed")
        call check(abs(norm_value - expected_L2) < TOL, "L2 norm: should equal 0.5")
    end subroutine test_density_norm_L2

    !> Test L∞ norm calculation
    !!
    !! Physics: L∞ norm = max|Δn(i)| finds the site with largest density change.
    !! Useful for detecting local convergence issues (e.g., edge states).
    subroutine test_density_norm_Linf()
        use fortuno_serial, only: check => serial_check
        use convergence_monitor, only: compute_density_norm, Linf
        integer, parameter :: L = 5
        real(dp) :: delta_n(L), norm_value
        real(dp) :: expected_Linf
        integer :: ierr

        ! Example: Δn = [0.1, -0.5, 0.2, 0.3, -0.4]
        delta_n = [0.1_dp, -0.5_dp, 0.2_dp, 0.3_dp, -0.4_dp]
        expected_Linf = 0.5_dp  ! max is at site 2

        call compute_density_norm(delta_n, L, Linf, norm_value, ierr)

        call check(ierr == 0, "Linf norm: computation should succeed")
        call check(abs(norm_value - expected_Linf) < TOL, &
                   "Linf norm: should equal 0.5 (max absolute value)")
    end subroutine test_density_norm_Linf

    !> Test invalid norm type handling
    !!
    !! Physics: Only L1, L2, Linf are defined. Other values are errors.
    subroutine test_density_norm_invalid_type()
        use fortuno_serial, only: check => serial_check
        use convergence_monitor, only: compute_density_norm
        use lsda_errors, only: ERROR_INVALID_INPUT
        integer, parameter :: L = 3
        real(dp) :: delta_n(L), norm_value
        integer :: ierr

        delta_n = [1.0_dp, 2.0_dp, 3.0_dp]

        call compute_density_norm(delta_n, L, 99, norm_value, ierr)  ! Invalid type

        call check(ierr == ERROR_INVALID_INPUT, &
                   "Invalid norm type: should return ERROR_INVALID_INPUT")
    end subroutine test_density_norm_invalid_type

    !> Test convergence check when converged
    !!
    !! Physics: SCF converged when ||Δn||₂ < tolerance (typically 1e-6).
    !! This stops the iteration and signals density self-consistency.
    subroutine test_convergence_check_converged()
        use fortuno_serial, only: check => serial_check
        use convergence_monitor, only: check_scf_convergence
        integer, parameter :: L = 5
        real(dp) :: delta_n(L)
        logical :: is_converged
        integer :: ierr

        ! Small differences: converged
        delta_n = [1.0e-7_dp, -2.0e-7_dp, 5.0e-8_dp, 1.0e-8_dp, -3.0e-8_dp]

        call check_scf_convergence(delta_n, L, 1.0e-6_dp, is_converged, ierr)

        call check(ierr == 0, "Converged: check should succeed")
        call check(is_converged, "Converged: small ||Δn|| should indicate convergence")
    end subroutine test_convergence_check_converged

    !> Test convergence check when not converged
    !!
    !! Physics: SCF not converged when ||Δn||₂ ≥ tolerance.
    !! Must continue iterating to achieve self-consistency.
    subroutine test_convergence_check_not_converged()
        use fortuno_serial, only: check => serial_check
        use convergence_monitor, only: check_scf_convergence
        integer, parameter :: L = 4
        real(dp) :: delta_n(L)
        logical :: is_converged
        integer :: ierr

        ! Large differences: not converged
        delta_n = [0.01_dp, -0.02_dp, 0.015_dp, -0.005_dp]

        call check_scf_convergence(delta_n, L, 1.0e-6_dp, is_converged, ierr)

        call check(ierr == 0, "Not converged: check should succeed")
        call check(.not. is_converged, &
                   "Not converged: large ||Δn|| should indicate non-convergence")
    end subroutine test_convergence_check_not_converged

    !> Test convergence check with custom tolerance
    !!
    !! Physics: Different problems need different tolerances.
    !! Tight tolerance (1e-8) for accurate energies, loose (1e-4) for quick tests.
    subroutine test_convergence_check_custom_tolerance()
        use fortuno_serial, only: check => serial_check
        use convergence_monitor, only: check_scf_convergence
        integer, parameter :: L = 3
        real(dp) :: delta_n(L)
        logical :: is_converged_tight, is_converged_loose
        integer :: ierr

        ! Medium-sized difference
        delta_n = [1.0e-5_dp, -1.0e-5_dp, 0.5e-5_dp]

        ! Tight tolerance: should NOT converge
        call check_scf_convergence(delta_n, L, 1.0e-6_dp, is_converged_tight, ierr)
        call check(.not. is_converged_tight, &
                   "Custom tol: tight tolerance should not converge")

        ! Loose tolerance: should converge
        call check_scf_convergence(delta_n, L, 1.0e-4_dp, is_converged_loose, ierr)
        call check(is_converged_loose, &
                   "Custom tol: loose tolerance should converge")
    end subroutine test_convergence_check_custom_tolerance

    !> Test convergence history initialization and cleanup
    !!
    !! Physics: Tracking convergence history allows analysis of SCF behavior:
    !! oscillations, monotonic convergence, sudden jumps (phase transitions).
    subroutine test_history_init_cleanup()
        use fortuno_serial, only: check => serial_check
        use convergence_monitor, only: convergence_history_t, &
                                       init_convergence_history, &
                                       cleanup_convergence_history
        type(convergence_history_t) :: history
        integer :: ierr

        ! Initialize for 10 iterations
        call init_convergence_history(history, 10, ierr)

        call check(ierr == 0, "History init: should succeed")
        call check(history%max_iter == 10, "History init: max_iter should be 10")
        call check(history%current_iter == 0, "History init: current_iter should be 0")
        call check(allocated(history%density_norms), &
                   "History init: density_norms should be allocated")
        call check(allocated(history%energies), &
                   "History init: energies should be allocated")

        ! Cleanup
        call cleanup_convergence_history(history, ierr)

        call check(ierr == 0, "History cleanup: should succeed")
        call check(.not. allocated(history%density_norms), &
                   "History cleanup: density_norms should be deallocated")
        call check(.not. allocated(history%energies), &
                   "History cleanup: energies should be deallocated")
    end subroutine test_history_init_cleanup

    !> Test updating convergence history
    !!
    !! Physics: Each SCF iteration produces a density norm and energy.
    !! Storing these allows monitoring convergence rate and detecting problems.
    subroutine test_history_update()
        use fortuno_serial, only: check => serial_check
        use convergence_monitor, only: convergence_history_t, &
                                       init_convergence_history, &
                                       update_convergence_history, &
                                       cleanup_convergence_history
        type(convergence_history_t) :: history
        integer :: ierr

        call init_convergence_history(history, 5, ierr)

        ! Update iteration 1
        call update_convergence_history(1, 1.0e-2_dp, -10.5_dp, history, ierr)
        call check(ierr == 0, "History update 1: should succeed")
        call check(history%current_iter == 1, "History update 1: current_iter = 1")
        call check(abs(history%density_norms(1) - 1.0e-2_dp) < TOL, &
                   "History update 1: norm stored correctly")
        call check(abs(history%energies(1) - (-10.5_dp)) < TOL, &
                   "History update 1: energy stored correctly")

        ! Update iteration 2
        call update_convergence_history(2, 5.0e-3_dp, -10.6_dp, history, ierr)
        call check(ierr == 0, "History update 2: should succeed")
        call check(history%current_iter == 2, "History update 2: current_iter = 2")

        call cleanup_convergence_history(history, ierr)
    end subroutine test_history_update

    !> Test bounds checking in history update
    !!
    !! Physics: Iteration number must be valid (1 ≤ iter ≤ max_iter).
    !! Out-of-bounds indicates programming error or corrupted state.
    subroutine test_history_bounds_checking()
        use fortuno_serial, only: check => serial_check
        use convergence_monitor, only: convergence_history_t, &
                                       init_convergence_history, &
                                       update_convergence_history, &
                                       cleanup_convergence_history
        use lsda_errors, only: ERROR_INVALID_INPUT
        type(convergence_history_t) :: history
        integer :: ierr

        call init_convergence_history(history, 3, ierr)

        ! Try iteration 0 (invalid)
        call update_convergence_history(0, 1.0_dp, -1.0_dp, history, ierr)
        call check(ierr == ERROR_INVALID_INPUT, &
                   "Bounds: iteration 0 should return ERROR_INVALID_INPUT")

        ! Try iteration 4 (beyond max_iter=3)
        call update_convergence_history(4, 1.0_dp, -1.0_dp, history, ierr)
        call check(ierr == ERROR_INVALID_INPUT, &
                   "Bounds: iteration > max_iter should return ERROR_INVALID_INPUT")

        call cleanup_convergence_history(history, ierr)
    end subroutine test_history_bounds_checking

end program test_convergence_monitor
