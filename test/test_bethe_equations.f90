program test_bethe_equations
    use fortuno_serial, only: execute_serial_cmd_app
    implicit none
    
    call execute_serial_cmd_app(get_bethe_tests())
    
contains

    function get_bethe_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests
        
        tests = test_list([ &
            test("theta_at_zero", test_theta_at_zero), &
            test("theta_antisymmetry", test_theta_antisymmetry), &
            test("Theta_at_zero", test_Theta_at_zero), &
            test("Theta_antisymmetry", test_Theta_antisymmetry), &
            test("dtheta_dx_numerical", test_dtheta_dx_numerical), &
            test("dTheta_capital_dx_numerical", test_dTheta_dx_numerical), &
            test("quantum_numbers_odd", test_quantum_numbers_odd), &
            test("quantum_numbers_even", test_quantum_numbers_even), &
            test("residual_dimensions", test_residual_dimensions), &
            test("jacobian_dimensions", test_jacobian_dimensions), &
            test("jacobian_diagonal", test_jacobian_diagonal) &
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

    subroutine test_theta_antisymmetry()
        use fortuno_serial, only: check => serial_check
        use bethe_equations, only: theta
        use lsda_constants, only: dp
        real(dp) :: x, U, theta_pos, theta_neg
        
        U = 4.0_dp
        x = 1.5_dp
        
        theta_pos = theta(x, U)
        theta_neg = theta(-x, U)
        
        call check(abs(theta_pos + theta_neg) < 1.0e-14_dp)
    end subroutine

    subroutine test_Theta_capital_at_zero()
        use fortuno_serial, only: check => serial_check
        use bethe_equations, only: Theta_capital
        use lsda_constants, only: dp
        real(dp) :: result, expected, U
        
        U = 4.0_dp
        result = Theta_capital(0.0_dp, U)
        expected = 0.0_dp
        
        call check(abs(result - expected) < 1.0e-14_dp)
    end subroutine

    subroutine test_Theta_capital_antisymmetry()
        use fortuno_serial, only: check => serial_check
        use bethe_equations, only: Theta_capital
        use lsda_constants, only: dp
        real(dp) :: x, U, Theta_pos, Theta_neg
        
        U = 4.0_dp
        x = 1.5_dp
        
        Theta_pos = Theta_capital(x, U)
        Theta_neg = Theta_capital(-x, U)
        
        call check(abs(Theta_pos + Theta_neg) < 1.0e-14_dp)
    end subroutine

    subroutine test_dtheta_dx_numerical()
        use fortuno_serial, only: check => serial_check
        use bethe_equations, only: theta, dtheta_dx
        use lsda_constants, only: dp
        real(dp) :: x, U, h, analytical, numerical
        
        U = 4.0_dp
        x = 1.5_dp
        h = 1.0e-6_dp

        ! Analytical derivative
        analytical = dtheta_dx(x, U)

        ! Numerical derivative (central finite differences)
        numerical = (theta(x + h, U) - theta(x - h, U)) / (2.0_dp * h)
        
        call check(abs(analytical - numerical) < 1.0e-8_dp)
    end subroutine

    subroutine test_dTheta_capital_dx_numerical()
        use fortuno_serial, only: check => serial_check
        use bethe_equations, only: Theta_capital, dTheta_capital_dx
        use lsda_constants, only: dp
        real(dp) :: x, U, h, analytical, numerical
        
        U = 4.0_dp
        x = 1.5_dp
        h = 1.0e-6_dp

        ! Analytical derivative
        analytical = dTheta_capital_dx(x, U)

        ! Numerical derivative (central finite differences)
        numerical = (Theta_capital(x + h, U) - Theta_capital(x - h, U)) / (2.0_dp * h)
        
        call check(abs(analytical - numerical) < 1.0e-8_dp)
    end subroutine

    subroutine test_quantum_numbers_odd()
        use fortuno_serial, only: check => serial_check
        use bethe_equations, only: initialize_quantum_numbers
        use lsda_constants, only: dp
        integer :: Nup, M
        real(dp), allocatable :: I(:), J(:)
        
        Nup = 5
        M = 3
        allocate(I(Nup), J(M))
        
        call initialize_quantum_numbers(Nup, M, I, J)

        ! Verify if I = [-2, -1, 0, 1, 2]
        call check(abs(I(1) - (-2.0_dp)) < 1.0e-14_dp)
        call check(abs(I(3) - 0.0_dp) < 1.0e-14_dp)
        call check(abs(I(5) - 2.0_dp) < 1.0e-14_dp)

        ! Verify if J = [-1, 0, 1]
        call check(abs(J(1) - (-1.0_dp)) < 1.0e-14_dp)
        call check(abs(J(2) - 0.0_dp) < 1.0e-14_dp)
        call check(abs(J(3) - 1.0_dp) < 1.0e-14_dp)
        
        deallocate(I, J)
    end subroutine

    subroutine test_quantum_numbers_even()
        use fortuno_serial, only: check => serial_check
        use bethe_equations, only: initialize_quantum_numbers
        use lsda_constants, only: dp
        integer :: Nup, M
        real(dp), allocatable :: I(:), J(:)
        
        Nup = 4
        M = 2
        allocate(I(Nup), J(M))
        
        call initialize_quantum_numbers(Nup, M, I, J)

        ! Verify if I = [-1.5, -0.5, 0.5, 1.5]
        call check(abs(I(1) - (-1.5_dp)) < 1.0e-14_dp)
        call check(abs(I(2) - (-0.5_dp)) < 1.0e-14_dp)
        call check(abs(I(3) - 0.5_dp) < 1.0e-14_dp)
        call check(abs(I(4) - 1.5_dp) < 1.0e-14_dp)
        
        deallocate(I, J)
    end subroutine

    subroutine test_residual_dimensions()
        use fortuno_serial, only: check => serial_check
        use bethe_equations, only: initialize_quantum_numbers, compute_residual
        use lsda_constants, only: dp, TWOPI
        integer :: Nup, M, L
        real(dp) :: U
        real(dp), allocatable :: k(:), Lambda(:), I(:), J(:), F(:)
        
        Nup = 3
        M = 2
        L = 10
        U = 4.0_dp
        
        allocate(k(Nup), Lambda(M), I(Nup), J(M))
        
        call initialize_quantum_numbers(Nup, M, I, J)
        
        ! Simple initial guess
        k = TWOPI * I / real(L, dp)
        Lambda = 0.0_dp

        ! Compute residual
        F = compute_residual(k, Lambda, I, J, L, U)

        ! Check correct dimension
        call check(size(F) == Nup + M)
        
        deallocate(k, Lambda, I, J)
    end subroutine

    subroutine test_jacobian_dimensions()
        use fortuno_serial, only: check => serial_check
        use bethe_equations, only: compute_jacobian
        use lsda_constants, only: dp
        integer :: Nup, M, L, total_size
        real(dp) :: U
        real(dp), allocatable :: k(:), Lambda(:), J_matrix(:,:)
        
        Nup = 3
        M = 2
        L = 10
        U = 4.0_dp
        total_size = Nup + M
        
        allocate(k(Nup), Lambda(M))
        
        k = [0.5_dp, 1.0_dp, -0.5_dp]
        Lambda = [0.2_dp, -0.3_dp]
        
        J_matrix = compute_jacobian(k, Lambda, L, U)

        ! Verify dimensions
        call check(size(J_matrix, 1) == total_size)
        call check(size(J_matrix, 2) == total_size)
        
        deallocate(k, Lambda)
    end subroutine

    subroutine test_jacobian_diagonal()
        use fortuno_serial, only: check => serial_check
        use bethe_equations, only: compute_jacobian
        use lsda_constants, only: dp
        integer :: Nup, M, L, i
        real(dp) :: U
        real(dp), allocatable :: k(:), Lambda(:), J_matrix(:,:)
        
        Nup = 3
        M = 2
        L = 10
        U = 4.0_dp
        
        allocate(k(Nup), Lambda(M))
        
        k = [0.5_dp, 1.0_dp, -0.5_dp]
        Lambda = [0.2_dp, -0.3_dp]
        
        J_matrix = compute_jacobian(k, Lambda, L, U)

        ! Verify that block A has approximately diagonal structure
        ! (diagonal dominant, since it has the term 1 - (1/L)*sum)
        ! Let's check that the diagonal is not zero
        do i = 1, Nup
            call check(abs(J_matrix(i, i)) > 0.5_dp)  ! Should be close to 1
        end do
        
        deallocate(k, Lambda)
    end subroutine

end program