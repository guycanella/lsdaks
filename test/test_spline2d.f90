!> Unit tests for spline2d module
program test_spline2d
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp
    implicit none

    real(dp), parameter :: TOL = 1.0e-9_dp

    call execute_serial_cmd_app(get_spline2d_tests())

contains

    function get_spline2d_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests

        tests = test_list([ &
            test("spline2d_init_destroy", test_spline2d_init_destroy), &
            test("spline2d_exact_at_grid", test_spline2d_exact_at_grid), &
            test("spline2d_linear_function", test_spline2d_linear_function), &
            test("spline2d_separable_function", test_spline2d_separable_function), &
            test("spline2d_interpolation_bounds", test_spline2d_interpolation_bounds) &
        ])
    end function get_spline2d_tests

    subroutine test_spline2d_linear_function()
        use fortuno_serial, only: check => serial_check
        use spline2d, only: spline2d_t, spline2d_init, spline2d_eval, spline2d_destroy
        use lsda_constants, only: dp

        type(spline2d_t) :: spl
        integer, parameter :: nx = 5, ny_max = 4
        real(dp) :: x_grid(nx), y_grid(ny_max, nx), f_values(ny_max, nx)
        integer :: n_y_pts(nx)
        real(dp) :: x_test, y_test, f_interp, f_exact, error
        integer :: i, j

        do i = 1, nx
            x_grid(i) = real(i, dp) * 0.2_dp
            n_y_pts(i) = ny_max
            do j = 1, ny_max
                y_grid(j, i) = real(j, dp) * 0.15_dp
                f_values(j, i) = 2.0_dp * x_grid(i) + 3.0_dp * y_grid(j, i) + 1.0_dp
            end do
        end do

        call spline2d_init(spl, x_grid, y_grid, f_values, n_y_pts)

        x_test = 0.45_dp
        y_test = 0.32_dp
        f_interp = spline2d_eval(spl, x_test, y_test)
        f_exact = 2.0_dp * x_test + 3.0_dp * y_test + 1.0_dp
        error = abs(f_interp - f_exact)

        call check(error < 1.0e-8_dp, "2D spline should be exact for linear functions")

        call spline2d_destroy(spl)

    end subroutine test_spline2d_linear_function

    subroutine test_spline2d_interpolation_bounds()
        use fortuno_serial, only: check => serial_check
        use spline2d, only: spline2d_t, spline2d_init, spline2d_eval, spline2d_destroy
        use lsda_constants, only: dp

        type(spline2d_t) :: spl
        integer, parameter :: nx = 1, ny_max = 5
        real(dp) :: x_grid(nx), y_grid(ny_max, nx), f_values(ny_max, nx)
        integer :: n_y_pts(nx)
        real(dp) :: y_test, f_interp
        integer :: j

        x_grid(1) = 0.5_dp
        n_y_pts(1) = ny_max

        do j = 1, ny_max
            y_grid(j, 1) = real(j-1, dp) * 0.25_dp
            f_values(j, 1) = sin(y_grid(j, 1))
        end do

        call spline2d_init(spl, x_grid, y_grid, f_values, n_y_pts)

        y_test = 0.35_dp
        f_interp = spline2d_eval(spl, x_grid(1), y_test)

        call check(f_interp == f_interp, "Should return valid number for single x point")

        call spline2d_destroy(spl)

    end subroutine test_spline2d_interpolation_bounds

    subroutine test_spline2d_init_destroy()
        use fortuno_serial, only: check => serial_check
        use spline2d, only: spline2d_t, spline2d_init, spline2d_destroy
        use lsda_constants, only: dp

        type(spline2d_t) :: spl
        integer, parameter :: nx = 3, ny_max = 4
        real(dp) :: x_grid(nx), y_grid(ny_max, nx), f_values(ny_max, nx)
        integer :: n_y_pts(nx)
        integer :: i, j

        do i = 1, nx
            x_grid(i) = real(i, dp)
            n_y_pts(i) = ny_max
            do j = 1, ny_max
                y_grid(j, i) = real(j, dp) * 0.5_dp
                f_values(j, i) = x_grid(i) + y_grid(j, i)
            end do
        end do

        call spline2d_init(spl, x_grid, y_grid, f_values, n_y_pts)

        call check(spl%initialized, "Spline should be initialized")
        call check(spl%n_x == nx, "Should have correct n_x")
        call check(spl%n_y_max == ny_max, "Should have correct n_y_max")
        call check(allocated(spl%x), "x array should be allocated")
        call check(allocated(spl%d2f_dy2), "d2f_dy2 array should be allocated")

        call spline2d_destroy(spl)
        call check(.not. spl%initialized, "Spline should be deinitialized")
        call check(.not. allocated(spl%x), "x array should be deallocated")

    end subroutine test_spline2d_init_destroy

    subroutine test_spline2d_exact_at_grid()
        use fortuno_serial, only: check => serial_check
        use spline2d, only: spline2d_t, spline2d_init, spline2d_eval, spline2d_destroy
        use lsda_constants, only: dp

        type(spline2d_t) :: spl
        integer, parameter :: nx = 4, ny_max = 5
        real(dp) :: x_grid(nx), y_grid(ny_max, nx), f_values(ny_max, nx)
        integer :: n_y_pts(nx)
        real(dp) :: f_interp, error
        integer :: i, j

        do i = 1, nx
            x_grid(i) = real(i, dp) * 0.3_dp
            n_y_pts(i) = ny_max
            do j = 1, ny_max
                y_grid(j, i) = real(j, dp) * 0.2_dp
                f_values(j, i) = x_grid(i) + 2.0_dp * y_grid(j, i)
            end do
        end do

        call spline2d_init(spl, x_grid, y_grid, f_values, n_y_pts)

        do i = 1, nx
            do j = 1, n_y_pts(i)
                f_interp = spline2d_eval(spl, x_grid(i), y_grid(j, i))
                error = abs(f_interp - f_values(j, i))
                call check(error < TOL, "2D spline should be exact at grid points")
            end do
        end do

        call spline2d_destroy(spl)

    end subroutine test_spline2d_exact_at_grid

    subroutine test_spline2d_separable_function()
        use fortuno_serial, only: check => serial_check
        use spline2d, only: spline2d_t, spline2d_init, spline2d_eval, spline2d_destroy
        use lsda_constants, only: dp

        type(spline2d_t) :: spl
        integer, parameter :: nx = 5, ny_max = 6
        real(dp) :: x_grid(nx), y_grid(ny_max, nx), f_values(ny_max, nx)
        integer :: n_y_pts(nx)
        real(dp) :: x_test, y_test, f_interp, f_exact, error
        integer :: i, j

        do i = 1, nx
            x_grid(i) = real(i-1, dp) * 0.25_dp
            n_y_pts(i) = ny_max
            do j = 1, ny_max
                y_grid(j, i) = real(j-1, dp) * 0.2_dp
                f_values(j, i) = (x_grid(i) + 1.0_dp) * (y_grid(j, i) + 2.0_dp)
            end do
        end do

        call spline2d_init(spl, x_grid, y_grid, f_values, n_y_pts)

        x_test = 0.35_dp
        y_test = 0.45_dp
        f_interp = spline2d_eval(spl, x_test, y_test)
        f_exact = (x_test + 1.0_dp) * (y_test + 2.0_dp)
        error = abs(f_interp - f_exact)

        call check(error < 1.0e-4_dp, "2D spline should accurately interpolate separable functions")

        call spline2d_destroy(spl)

    end subroutine test_spline2d_separable_function

end program test_spline2d
