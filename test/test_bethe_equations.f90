program test_bethe_equations
    use fortuno_serial, only: execute_serial_cmd_app
    implicit none
    
    call execute_serial_cmd_app(get_bethe_tests())
    
contains

    function get_bethe_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests
        
        tests = test_list([ &
            test("theta_at_zero", test_theta_at_zero) &
        ])
    end function

    subroutine test_theta_at_zero()
        use fortuno_serial, only: check => serial_check
        use bethe_equations, only: theta
        use lsda_constants, only: dp
        real(dp) :: result, expected, U
        
        U = 4.0_dp
        result = theta(0.0_dp, U)
        expected = 0.0_dp
        
        call check(abs(result - expected) < 1.0e-14_dp)
    end subroutine

end program