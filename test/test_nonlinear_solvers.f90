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
            test("solve_inputs_not_modified", test_inputs_not_modified) &
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

end program