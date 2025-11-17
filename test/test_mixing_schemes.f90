!> Unit tests for mixing_schemes module
program test_mixing_schemes
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp
    implicit none

    real(dp), parameter :: TOL = 1.0e-9_dp

    call execute_serial_cmd_app(get_mixing_tests())

contains

    function get_mixing_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests

        tests = test_list([ &
            test("linear_mixing_alpha_half", test_linear_mixing_alpha_half), &
            test("linear_mixing_alpha_one", test_linear_mixing_alpha_one), &
            test("linear_mixing_alpha_small", test_linear_mixing_alpha_small), &
            test("linear_mixing_bounds", test_linear_mixing_bounds), &
            test("linear_mixing_convergence", test_linear_mixing_convergence), &
            test("linear_mixing_invalid_alpha_zero", test_linear_mixing_invalid_alpha_zero), &
            test("linear_mixing_invalid_alpha_negative", test_linear_mixing_invalid_alpha_negative), &
            test("linear_mixing_invalid_alpha_large", test_linear_mixing_invalid_alpha_large), &
            test("linear_mixing_size_mismatch", test_linear_mixing_size_mismatch) &
        ])
    end function get_mixing_tests

    !> Test linear mixing with α = 0.5 (equal weights)
    !!
    !! Physics: α = 0.5 gives n_mixed = 0.5·n_old + 0.5·n_new (simple average).
    !! Moderate damping, good starting point for most systems.
    subroutine test_linear_mixing_alpha_half()
        use fortuno_serial, only: check => serial_check
        use mixing_schemes, only: linear_mixing
        integer, parameter :: L = 4
        real(dp) :: n_new(L), n_old(L), n_mixed(L), expected(L)
        integer :: ierr, i

        n_old = [1.0_dp, 0.8_dp, 1.2_dp, 0.5_dp]
        n_new = [1.2_dp, 0.6_dp, 1.4_dp, 0.7_dp]
        expected = 0.5_dp * n_old + 0.5_dp * n_new

        call linear_mixing(n_new, n_old, 0.5_dp, n_mixed, L, ierr)

        call check(ierr == 0, "α=0.5: mixing should succeed")
        do i = 1, L
            call check(abs(n_mixed(i) - expected(i)) < TOL, &
                       "α=0.5: mixed density should be average")
        end do
    end subroutine test_linear_mixing_alpha_half

    !> Test linear mixing with α = 1.0 (full update)
    !!
    !! Physics: α = 1.0 gives n_mixed = n_new (no damping).
    !! Fastest convergence if SCF is stable, but can oscillate for strong correlations.
    subroutine test_linear_mixing_alpha_one()
        use fortuno_serial, only: check => serial_check
        use mixing_schemes, only: linear_mixing
        integer, parameter :: L = 3
        real(dp) :: n_new(L), n_old(L), n_mixed(L)
        integer :: ierr, i

        n_old = [0.5_dp, 1.0_dp, 1.5_dp]
        n_new = [0.7_dp, 1.2_dp, 1.3_dp]

        call linear_mixing(n_new, n_old, 1.0_dp, n_mixed, L, ierr)

        call check(ierr == 0, "α=1.0: mixing should succeed")
        do i = 1, L
            call check(abs(n_mixed(i) - n_new(i)) < TOL, &
                       "α=1.0: mixed density should equal n_new (no damping)")
        end do
    end subroutine test_linear_mixing_alpha_one

    !> Test linear mixing with small α = 0.1 (heavy damping)
    !!
    !! Physics: α = 0.1 gives n_mixed = 0.9·n_old + 0.1·n_new (slow update).
    !! Heavy damping prevents oscillations in difficult cases (e.g., Mott insulators).
    subroutine test_linear_mixing_alpha_small()
        use fortuno_serial, only: check => serial_check
        use mixing_schemes, only: linear_mixing
        integer, parameter :: L = 5
        real(dp) :: n_new(L), n_old(L), n_mixed(L), expected(L)
        integer :: ierr, i

        n_old = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
        n_new = [1.5_dp, 0.5_dp, 1.2_dp, 0.8_dp, 1.1_dp]
        expected = 0.9_dp * n_old + 0.1_dp * n_new

        call linear_mixing(n_new, n_old, 0.1_dp, n_mixed, L, ierr)

        call check(ierr == 0, "α=0.1: mixing should succeed")
        do i = 1, L
            call check(abs(n_mixed(i) - expected(i)) < TOL, &
                       "α=0.1: mixed density should be heavily damped")
        end do

        ! Check that n_mixed is much closer to n_old than to n_new
        call check(abs(n_mixed(1) - n_old(1)) < abs(n_mixed(1) - n_new(1)), &
                   "α=0.1: mixed should be closer to n_old (heavy damping)")
    end subroutine test_linear_mixing_alpha_small

    !> Test that mixing preserves physical density bounds
    !!
    !! Physics: If 0 ≤ n_old, n_new ≤ 2, then linear mixing guarantees
    !! 0 ≤ n_mixed ≤ 2 (convex combination preserves bounds).
    subroutine test_linear_mixing_bounds()
        use fortuno_serial, only: check => serial_check
        use mixing_schemes, only: linear_mixing
        integer, parameter :: L = 6
        real(dp) :: n_new(L), n_old(L), n_mixed(L)
        integer :: ierr, i

        ! Physical densities: 0 ≤ n ≤ 2
        n_old = [0.0_dp, 0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp, 1.2_dp]
        n_new = [0.2_dp, 0.8_dp, 0.9_dp, 1.8_dp, 1.9_dp, 0.7_dp]

        call linear_mixing(n_new, n_old, 0.3_dp, n_mixed, L, ierr)

        call check(ierr == 0, "Bounds: mixing should succeed")
        do i = 1, L
            call check(n_mixed(i) >= 0.0_dp - TOL, &
                       "Bounds: mixed density should be non-negative")
            call check(n_mixed(i) <= 2.0_dp + TOL, &
                       "Bounds: mixed density should be ≤ 2 (Pauli exclusion)")
        end do
    end subroutine test_linear_mixing_bounds

    !> Test convergence behavior with repeated mixing
    !!
    !! Physics: In SCF, if n_new converges to n* (fixed point), then
    !! repeated mixing n_mixed → n_new → n_mixed should approach n*.
    !! This tests the fundamental SCF convergence mechanism.
    subroutine test_linear_mixing_convergence()
        use fortuno_serial, only: check => serial_check
        use mixing_schemes, only: linear_mixing
        integer, parameter :: L = 3
        real(dp) :: n_old(L), n_new(L), n_mixed(L)
        real(dp) :: alpha = 0.5_dp
        integer :: ierr, iter

        ! Start with different densities
        n_old = [0.5_dp, 1.0_dp, 1.5_dp]
        n_new = [1.5_dp, 1.0_dp, 0.5_dp]

        ! Simulate SCF: assume n_new converges to fixed point [1.0, 1.0, 1.0]
        do iter = 1, 10
            call linear_mixing(n_new, n_old, alpha, n_mixed, L, ierr)
            n_old = n_mixed
            ! In real SCF, we'd diagonalize H(n_mixed) to get new n_new
            ! Here we simulate convergence to [1.0, 1.0, 1.0]
            n_new = n_new * 0.8_dp + [1.0_dp, 1.0_dp, 1.0_dp] * 0.2_dp
        end do

        ! After 10 iterations, should be closer to fixed point
        call check(abs(n_mixed(1) - 1.0_dp) < 0.5_dp, &
                   "Convergence: density should approach fixed point")
        call check(abs(n_mixed(2) - 1.0_dp) < 0.1_dp, &
                   "Convergence: site 2 already at fixed point")
    end subroutine test_linear_mixing_convergence

    !> Test invalid α = 0 (no update)
    !!
    !! Physics: α = 0 would give n_mixed = n_old (frozen density, no SCF progress).
    !! This is invalid as it prevents convergence.
    subroutine test_linear_mixing_invalid_alpha_zero()
        use fortuno_serial, only: check => serial_check
        use mixing_schemes, only: linear_mixing
        use lsda_errors, only: ERROR_INVALID_INPUT
        integer, parameter :: L = 3
        real(dp) :: n_new(L), n_old(L), n_mixed(L)
        integer :: ierr

        n_old = [1.0_dp, 1.0_dp, 1.0_dp]
        n_new = [1.1_dp, 0.9_dp, 1.0_dp]

        call linear_mixing(n_new, n_old, 0.0_dp, n_mixed, L, ierr)

        call check(ierr == ERROR_INVALID_INPUT, &
                   "Invalid α=0: should return ERROR_INVALID_INPUT")
    end subroutine test_linear_mixing_invalid_alpha_zero

    !> Test invalid α < 0 (unphysical)
    !!
    !! Physics: Negative α can lead to unphysical densities (negative or > 2).
    !! Must be rejected.
    subroutine test_linear_mixing_invalid_alpha_negative()
        use fortuno_serial, only: check => serial_check
        use mixing_schemes, only: linear_mixing
        use lsda_errors, only: ERROR_INVALID_INPUT
        integer, parameter :: L = 3
        real(dp) :: n_new(L), n_old(L), n_mixed(L)
        integer :: ierr

        n_old = [1.0_dp, 1.0_dp, 1.0_dp]
        n_new = [1.1_dp, 0.9_dp, 1.0_dp]

        call linear_mixing(n_new, n_old, -0.5_dp, n_mixed, L, ierr)

        call check(ierr == ERROR_INVALID_INPUT, &
                   "Invalid α<0: should return ERROR_INVALID_INPUT")
    end subroutine test_linear_mixing_invalid_alpha_negative

    !> Test invalid α > 1 (over-relaxation)
    !!
    !! Physics: α > 1 is over-relaxation, can destabilize SCF.
    !! While mathematically valid, it's dangerous and should be rejected
    !! for safety (unless user explicitly wants advanced schemes).
    subroutine test_linear_mixing_invalid_alpha_large()
        use fortuno_serial, only: check => serial_check
        use mixing_schemes, only: linear_mixing
        use lsda_errors, only: ERROR_INVALID_INPUT
        integer, parameter :: L = 3
        real(dp) :: n_new(L), n_old(L), n_mixed(L)
        integer :: ierr

        n_old = [1.0_dp, 1.0_dp, 1.0_dp]
        n_new = [1.1_dp, 0.9_dp, 1.0_dp]

        call linear_mixing(n_new, n_old, 1.5_dp, n_mixed, L, ierr)

        call check(ierr == ERROR_INVALID_INPUT, &
                   "Invalid α>1: should return ERROR_INVALID_INPUT")
    end subroutine test_linear_mixing_invalid_alpha_large

    !> Test size mismatch detection
    !!
    !! Physics: All density arrays must have same size L.
    !! Mismatch indicates programming error or memory corruption.
    subroutine test_linear_mixing_size_mismatch()
        use fortuno_serial, only: check => serial_check
        use mixing_schemes, only: linear_mixing
        use lsda_errors, only: ERROR_SIZE_MISMATCH
        integer, parameter :: L = 4
        real(dp) :: n_new(L), n_old(L), n_mixed(L+1)  ! Wrong size!
        integer :: ierr

        n_old = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
        n_new = [1.1_dp, 0.9_dp, 1.0_dp, 1.05_dp]

        call linear_mixing(n_new, n_old, 0.5_dp, n_mixed, L, ierr)

        call check(ierr == ERROR_SIZE_MISMATCH, &
                   "Size mismatch: should return ERROR_SIZE_MISMATCH")
    end subroutine test_linear_mixing_size_mismatch

end program test_mixing_schemes
