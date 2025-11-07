program test_nonlinear_solvers
    use fortuno_serial, only: execute_serial_cmd_app
    implicit none

    call execute_serial_cmd_app(get_nonlinear_solver_tests())

contains

    function get_nonlinear_solver_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests
        
        tests = test_list([ &
            test("solve_2x2", test_solve_2x2), &
            test("solve_identity", test_solve_identity), &
            test("solve_inputs_not_modified", test_inputs_not_modified), &
            test("line_search_full_step", test_line_search_full_step), &
            test("line_search_positive_alpha", test_line_search_positive_alpha) &
        ])
    end function

    subroutine test_solve_2x2()
        use fortuno_serial, only: check => serial_check
        use nonlinear_solvers, only: solve_linear_system
        use lsda_constants, only: dp
        
        real(dp) :: A(2,2), b(2), x(2), x_expected(2)
        
        A = reshape([2.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [2, 2])
        b = [5.0_dp, 7.0_dp]
        x_expected = [1.6_dp, 1.8_dp]
        
        call solve_linear_system(A, b, x)
        
        call check(abs(x(1) - x_expected(1)) < 1.0e-12_dp)
        call check(abs(x(2) - x_expected(2)) < 1.0e-12_dp)
    end subroutine

    subroutine test_solve_identity()
        use fortuno_serial, only: check => serial_check
        use nonlinear_solvers, only: solve_linear_system
        use lsda_constants, only: dp

        real(dp) :: A(2,2), b(2), x(2)

        ! If A = I (identity), then x must be equal to b
        A = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
        b = [3.0_dp, 7.0_dp]
        
        call solve_linear_system(A, b, x)
        
        call check(abs(x(1) - b(1)) < 1.0e-14_dp)
        call check(abs(x(2) - b(2)) < 1.0e-14_dp)
    end subroutine

    subroutine test_inputs_not_modified()
        use fortuno_serial, only: check => serial_check
        use nonlinear_solvers, only: solve_linear_system
        use lsda_constants, only: dp

        real(dp) :: A(2,2), A_original(2,2)
        real(dp) :: b(2), b_original(2)
        real(dp) :: x(2)
        
        A = reshape([2.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [2, 2])
        b = [5.0_dp, 7.0_dp]
        
        A_original = A
        b_original = b
        
        call solve_linear_system(A, b, x)

        ! Check that A and b have not changed
        call check(all(abs(A - A_original) < 1.0e-14_dp))
        call check(all(abs(b - b_original) < 1.0e-14_dp))
    end subroutine

    subroutine test_line_search_full_step()
        use fortuno_serial, only: check => serial_check
        use nonlinear_solvers, only: line_search
        use bethe_equations, only: initialize_quantum_numbers, compute_residual
        use lsda_constants, only: dp, TWOPI
        
        integer :: Nup, M, L
        real(dp) :: U, alpha
        real(dp), allocatable :: x(:), dx(:), F_old(:), I(:), J(:)
        
        ! Small system
        Nup = 2
        M = 1
        L = 10
        U = 4.0_dp
        
        allocate(x(Nup+M), dx(Nup+M), I(Nup), J(M))
        
        call initialize_quantum_numbers(Nup, M, I, J)
        
        ! Initial position
        x = [0.1_dp, -0.1_dp, 0.05_dp]

        ! Direction that reduces the residual (simulated)
        dx = -0.01_dp * x  ! Small step in the opposite direction

        F_old = compute_residual(x(1:Nup), x(Nup+1:), I, J, L, U)
        alpha = line_search(x, dx, F_old, I, J, L, U)

        ! For small and good direction, should accept full step
        call check(alpha > 0.5_dp)  ! Alpha should be reasonably large

        deallocate(x, dx, I, J)
    end subroutine

    subroutine test_line_search_positive_alpha()
        use fortuno_serial, only: check => serial_check
        use nonlinear_solvers, only: line_search
        use bethe_equations, only: initialize_quantum_numbers, compute_residual
        use lsda_constants, only: dp
        
        integer :: Nup, M, L
        real(dp) :: U, alpha
        real(dp), allocatable :: x(:), dx(:), F_old(:), I(:), J(:)
        
        Nup = 2
        M = 1
        L = 10
        U = 4.0_dp
        
        allocate(x(Nup+M), dx(Nup+M), I(Nup), J(M))
        
        call initialize_quantum_numbers(Nup, M, I, J)
        
        x = [0.5_dp, -0.3_dp, 0.2_dp]
        dx = [0.1_dp, 0.1_dp, 0.1_dp]
        
        F_old = compute_residual(x(1:Nup), x(Nup+1:), I, J, L, U)
        
        alpha = line_search(x, dx, F_old, I, J, L, U)
        
        ! Alpha should be always positive and less than or equal to 1
        call check(alpha > 0.0_dp)
        call check(alpha <= 1.0_dp)
        
        deallocate(x, dx, I, J)
    end subroutine

end program