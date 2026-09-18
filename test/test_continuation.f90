program test_continuation
    use fortuno_serial, only: execute_serial_cmd_app
    implicit none
    
    call execute_serial_cmd_app(get_continuation_tests())
    
contains

    function get_continuation_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests
        
        tests = test_list([ &
            test("estimate_dxdU_simple", test_estimate_dxdU_simple), &
            test("sweep_forward_3pts", test_sweep_forward_3pts), &
            test("sweep_forward_all_converge", test_sweep_forward_all_converge), &
            test("sweep_backward_3pts", test_sweep_backward_3pts), &
            test("sweep_bidirectional_consistency", test_bidirectional_consistency), &
            test("merge_sweeps_masks_failures", test_merge_sweeps_masks_failures) &
        ])
    end function


    !> Regression test for the forward/backward merge.
    !!
    !! `store_point` writes NaN at a point a sweep failed to solve.  The merge
    !! used to be an unconditional average, which had two consequences:
    !!
    !! 1. a point the forward sweep solved perfectly came out as NaN just
    !!    because the backward sweep failed at the same `U`;
    !! 2. `maxval(abs(E_fwd - E_bwd))` became NaN, and `NaN > 1e-6` is
    !!    `.false.` in IEEE arithmetic, so the forward/backward consistency
    !!    warning - the only cross-check in the module - never fired again for
    !!    the rest of the sweep.  A genuine disagreement at another point was
    !!    therefore reported as consistent.
    !!
    !! Both halves are checked here on synthetic sweeps, because driving
    !! `solve_newton` into failure from a test is not practical.
    subroutine test_merge_sweeps_masks_failures()
        use fortuno_serial, only: check => serial_check
        use lsda_constants, only: dp
        use continuation, only: merge_sweeps
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan

        integer, parameter :: NP = 3, NX = 2
        real(dp) :: sol_fwd(NX, NP), sol_bwd(NX, NP), E_fwd(NP), E_bwd(NP)
        real(dp) :: solutions(NX, NP), energies(NP), max_diff, nan
        logical :: flags_fwd(NP), flags_bwd(NP), flags(NP), inconsistent

        nan = ieee_value(0.0_dp, ieee_quiet_nan)

        ! Point 1: both directions fine and in agreement.
        ! Point 2: the backward sweep failed; the forward value is good.
        ! Point 3: both fine but disagreeing well above the 1e-6 tolerance.
        sol_fwd(:, 1) = [1.0_dp, 2.0_dp]
        sol_bwd(:, 1) = [1.0_dp, 2.0_dp]
        E_fwd(1) = -3.0_dp
        E_bwd(1) = -3.0_dp

        sol_fwd(:, 2) = [5.0_dp, 6.0_dp]
        sol_bwd(:, 2) = nan
        E_fwd(2) = -7.0_dp
        E_bwd(2) = nan

        sol_fwd(:, 3) = [9.0_dp, 10.0_dp]
        sol_bwd(:, 3) = [9.0_dp, 10.0_dp]
        E_fwd(3) = -11.0_dp
        E_bwd(3) = -11.001_dp

        flags_fwd = [.true., .true., .true.]
        flags_bwd = [.true., .false., .true.]

        call merge_sweeps(sol_fwd, E_fwd, flags_fwd, sol_bwd, E_bwd, flags_bwd, &
                          solutions, energies, flags, max_diff, inconsistent)

        ! (a) the good forward point survives the backward failure
        call check(.not. ieee_is_nan(energies(2)), &
                   "a point solved by one direction must not be poisoned by the other")
        call check(abs(energies(2) + 7.0_dp) < 1.0e-14_dp, &
                   "the surviving direction must be copied through unchanged")
        call check(all(abs(solutions(:, 2) - sol_fwd(:, 2)) < 1.0e-14_dp), &
                   "the surviving solution vector must be copied through unchanged")

        ! (b) the inconsistency warning still fires despite the NaN point
        call check(inconsistent, &
                   "forward/backward disagreement must be reported even when " // &
                   "another point failed")
        call check(abs(max_diff - 1.0e-3_dp) < 1.0e-12_dp, &
                   "max_diff must be measured over the cross-checked points only")

        ! averaging and flags on the fully cross-checked points
        call check(abs(energies(1) + 3.0_dp) < 1.0e-14_dp, &
                   "cross-checked points must be averaged")
        call check(flags(1) .and. flags(3) .and. .not. flags(2), &
                   "only points accepted by both directions count as converged")

        ! (c) a NaN that escapes the mask must fire the warning **on its own**.
        ! Point 3 is brought back into agreement, so the NaN at point 2 is the
        ! only evidence left; and it sits *before* point 3 in the loop, which
        ! is what a non-sticky `max_diff` cannot survive (the NaN it stored is
        ! overwritten by the next finite difference, because `.not. (diff <=
        ! NaN)` is `.true.`).  Without the accumulated NaN flag this check
        ! reports max_diff = 0 and inconsistent = .false.
        flags_bwd(2) = .true.
        E_bwd(3) = -11.0_dp
        call merge_sweeps(sol_fwd, E_fwd, flags_fwd, sol_bwd, E_bwd, flags_bwd, &
                          solutions, energies, flags, max_diff, inconsistent)
        call check(inconsistent, &
                   "a NaN energy difference must fire the warning even when " // &
                   "every other point agrees")
        call check(ieee_is_nan(energies(2)), &
                   "a point averaged with a NaN must stay NaN, not be reported clean")
    end subroutine


    subroutine test_estimate_dxdU_simple()
        use fortuno_serial, only: check => serial_check
        use continuation, only: estimate_dxdU
        use lsda_constants, only: dp
        
        real(dp) :: x_old(3), x_current(3), dU, dxdU(3), expected(3)
        
        dU = 2.0_dp
        x_old = [1.0_dp, 2.0_dp, 3.0_dp]
        x_current = [2.0_dp, 4.0_dp, 6.0_dp]
        
        dxdU = estimate_dxdU(x_old, x_current, dU)
        
        expected = [0.5_dp, 1.0_dp, 1.5_dp]
        
        call check(all(abs(dxdU - expected) < 1.0e-14_dp), &
                   "dxdU should match expected derivative")
    end subroutine


    subroutine test_sweep_forward_3pts()
        use fortuno_serial, only: check => serial_check
        use lsda_constants, only: dp
        use bethe_equations, only: initialize_quantum_numbers
        use continuation, only: sweep_U_forward
        
        integer :: Nup, M, L, n_points
        real(dp), allocatable :: I(:), J(:), U_values(:)
        real(dp), allocatable :: solutions(:,:), energies(:)
        logical, allocatable :: converged_flags(:)
        
        Nup = 3
        M = 2
        L = 10
        n_points = 3
        
        allocate(I(Nup), J(M))
        allocate(U_values(n_points))
        allocate(solutions(Nup+M, n_points))
        allocate(energies(n_points))
        allocate(converged_flags(n_points))
        
        call initialize_quantum_numbers(Nup, M, I, J)
        U_values = [0.0_dp, 2.0_dp, 4.0_dp]
        
        call sweep_U_forward(I, J, L, U_values, solutions, energies, converged_flags)
        
        call check(all(converged_flags), "All points should converge")
        
        call check(size(solutions, 1) == Nup + M, "Solution dimension 1")
        call check(size(solutions, 2) == n_points, "Solution dimension 2")
        call check(size(energies) == n_points, "Energy dimension")
        
        deallocate(I, J, U_values, solutions, energies, converged_flags)
    end subroutine


    subroutine test_sweep_forward_all_converge()
        use fortuno_serial, only: check => serial_check
        use lsda_constants, only: dp
        use bethe_equations, only: initialize_quantum_numbers
        use continuation, only: sweep_U_forward
        
        integer :: Nup, M, L, n_points, k
        real(dp), allocatable :: I(:), J(:), U_values(:)
        real(dp), allocatable :: solutions(:,:), energies(:)
        logical, allocatable :: converged_flags(:)
        
        Nup = 2
        M = 1
        L = 10
        n_points = 5
        
        allocate(I(Nup), J(M))
        allocate(U_values(n_points))
        allocate(solutions(Nup+M, n_points))
        allocate(energies(n_points))
        allocate(converged_flags(n_points))
        
        call initialize_quantum_numbers(Nup, M, I, J)
        
        U_values = [(real(k, dp), k=0, n_points-1)]
        
        call sweep_U_forward(I, J, L, U_values, solutions, energies, converged_flags)
        
        call check(all(converged_flags), "All 5 points should converge")
        
        call check(all(energies < 0.0_dp), "All energies should be negative")
        
        deallocate(I, J, U_values, solutions, energies, converged_flags)
    end subroutine


    subroutine test_sweep_backward_3pts()
        use fortuno_serial, only: check => serial_check
        use lsda_constants, only: dp
        use bethe_equations, only: initialize_quantum_numbers
        use continuation, only: sweep_U_backward
        
        integer :: Nup, M, L, n_points
        real(dp), allocatable :: I(:), J(:), U_values(:)
        real(dp), allocatable :: solutions(:,:), energies(:)
        logical, allocatable :: converged_flags(:)
        
        Nup = 3
        M = 2
        L = 10
        n_points = 3
        
        allocate(I(Nup), J(M))
        allocate(U_values(n_points))
        allocate(solutions(Nup+M, n_points))
        allocate(energies(n_points))
        allocate(converged_flags(n_points))
        
        call initialize_quantum_numbers(Nup, M, I, J)
        U_values = [0.0_dp, 2.0_dp, 4.0_dp]
        
        call sweep_U_backward(I, J, L, U_values, solutions, energies, converged_flags)
        
        call check(all(converged_flags), "All points should converge (backward)")
        
        deallocate(I, J, U_values, solutions, energies, converged_flags)
    end subroutine


    subroutine test_bidirectional_consistency()
        use fortuno_serial, only: check => serial_check
        use lsda_constants, only: dp
        use bethe_equations, only: initialize_quantum_numbers
        use continuation, only: sweep_U_forward, sweep_U_backward, sweep_U_bidirectional
        
        integer :: Nup, M, L, n_points
        real(dp), allocatable :: I(:), J(:), U_values(:)
        real(dp), allocatable :: sol_fwd(:,:), sol_bwd(:,:), sol_bi(:,:)
        real(dp), allocatable :: E_fwd(:), E_bwd(:), E_bi(:)
        logical, allocatable :: flags_fwd(:), flags_bwd(:), flags_bi(:)
        real(dp) :: max_diff
        
        Nup = 2
        M = 1
        L = 10
        n_points = 5
        
        allocate(I(Nup), J(M))
        allocate(U_values(n_points))
        allocate(sol_fwd(Nup+M, n_points), sol_bwd(Nup+M, n_points), sol_bi(Nup+M, n_points))
        allocate(E_fwd(n_points), E_bwd(n_points), E_bi(n_points))
        allocate(flags_fwd(n_points), flags_bwd(n_points), flags_bi(n_points))
        
        call initialize_quantum_numbers(Nup, M, I, J)
        U_values = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
        
        call sweep_U_forward(I, J, L, U_values, sol_fwd, E_fwd, flags_fwd)
        
        call sweep_U_backward(I, J, L, U_values, sol_bwd, E_bwd, flags_bwd)
        
        call sweep_U_bidirectional(I, J, L, U_values, sol_bi, E_bi, flags_bi)
        
        max_diff = maxval(abs(E_fwd - E_bwd))
        call check(max_diff < 1.0e-4_dp, &
                   "Forward and backward energies should agree")
        
        call check(all(abs(E_bi - (E_fwd + E_bwd)/2.0_dp) < 1.0e-14_dp), &
                   "Bidirectional should be average of forward and backward")
        
        call check(all(flags_fwd), "Forward sweep should converge")
        call check(all(flags_bwd), "Backward sweep should converge")
        call check(all(flags_bi), "Bidirectional sweep should converge")
        
        deallocate(I, J, U_values)
        deallocate(sol_fwd, sol_bwd, sol_bi)
        deallocate(E_fwd, E_bwd, E_bi)
        deallocate(flags_fwd, flags_bwd, flags_bi)
    end subroutine

end program test_continuation