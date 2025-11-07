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
            test("sweep_bidirectional_consistency", test_bidirectional_consistency) &
        ])
    end function


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