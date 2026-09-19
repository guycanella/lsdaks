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
            test("spline2d_uninitialized_guard", test_spline2d_uninitialized_guard), &
            test("spline2d_exact_at_grid", test_spline2d_exact_at_grid), &
            test("spline2d_linear_function", test_spline2d_linear_function), &
            test("spline2d_separable_function", test_spline2d_separable_function), &
            test("spline2d_cubic_in_x", test_spline2d_cubic_in_x), &
            test("spline2d_synthetic_node", test_spline2d_synthetic_node), &
            test("spline2d_interpolation_bounds", test_spline2d_interpolation_bounds) &
        ])
    end function get_spline2d_tests

    !> An uninitialized spline must return before any dimension or data array is used.
    !!
    !! A default spline has zero metadata and no allocated data arrays. The
    !! guard must return before accessing its data; this status is the
    !! observable contract.
    subroutine test_spline2d_uninitialized_guard()
        use fortuno_serial, only: check => serial_check
        use spline2d, only: spline2d_t, spline2d_eval
        use lsda_errors, only: ERROR_NOT_INITIALIZED

        type(spline2d_t) :: spl
        real(dp) :: value
        integer :: ierr

        ! The default initialized state is .false.; this is deterministic
        ! valid input to the public error path.
        value = spline2d_eval(spl, 0.0_dp, 0.0_dp, ierr = ierr)
        call check(ierr == ERROR_NOT_INITIALIZED, &
                   "uninitialized spline must report ERROR_NOT_INITIALIZED")
        call check(abs(value) < TOL, "uninitialized spline should return its guarded value")
    end subroutine test_spline2d_uninitialized_guard

    !> The x direction must be interpolated with a cubic spline, not linearly
    !!
    !! Samples the separable function f(x, y) = sin(x)·(2 + y) on a coarse grid
    !! and checks the interpolation error at interior points that are not grid
    !! nodes. The y dependence is deliberately LINEAR: the y-direction spline
    !! reproduces it exactly, so the whole residual error comes from the x
    !! direction and the test measures precisely what it claims to.
    !!
    !! (The task sheet suggested sin(x)·y², but the y-direction spline uses
    !! natural boundary conditions, which alone cost ~1e-2 on a quadratic and
    !! would swamp the x-direction error under test.)
    !!
    !! With the previous linear interpolation in x the error at these points is
    !! ~1e-2, i.e. four orders of magnitude above the threshold below.
    subroutine test_spline2d_cubic_in_x()
        use fortuno_serial, only: check => serial_check
        use spline2d, only: spline2d_t, spline2d_init, spline2d_eval, spline2d_destroy
        use lsda_constants, only: dp

        type(spline2d_t) :: spl
        integer, parameter :: nx = 11, ny_max = 6
        real(dp) :: x_grid(nx), y_grid(ny_max, nx), f_values(ny_max, nx)
        integer :: n_y_pts(nx)
        real(dp) :: x_test, y_test, f_interp, f_exact, error
        integer :: i, j, k

        real(dp), parameter :: X_PROBE(4) = [0.5_dp, 0.7_dp, 0.9_dp, 1.1_dp]
        real(dp), parameter :: Y_PROBE(2) = [0.35_dp, 0.6_dp]

        do i = 1, nx
            x_grid(i) = real(i-1, dp) * 0.2_dp
            n_y_pts(i) = ny_max
            do j = 1, ny_max
                y_grid(j, i) = real(j-1, dp) * 0.25_dp
                f_values(j, i) = sin(x_grid(i)) * (2.0_dp + y_grid(j, i))
            end do
        end do

        call spline2d_init(spl, x_grid, y_grid, f_values, n_y_pts)

        do k = 1, size(X_PROBE)
            do j = 1, size(Y_PROBE)
                x_test = X_PROBE(k)
                y_test = Y_PROBE(j)
                f_interp = spline2d_eval(spl, x_test, y_test)
                f_exact = sin(x_test) * (2.0_dp + y_test)
                error = abs(f_interp - f_exact)
                call check(error < 1.0e-4_dp, &
                           "cubic interpolation in x should be accurate off the nodes")
            end do
        end do

        call spline2d_destroy(spl)

    end subroutine test_spline2d_cubic_in_x

    !> The synthetic node removes the extrapolation below the first abscissa
    !!
    !! When `node0_value` is supplied, a node is prepended at x = y and the row
    !! window starts at the last row whose abscissa is still >= y. Evaluating at
    !! x = y must therefore return the supplied value exactly, whatever the
    !! tabulated rows say, and no evaluation with x >= y can fall outside the
    !! local spline support.
    !!
    !! Without this mechanism a point with x below the first tabulated abscissa
    !! (the n = 0.0015 sites of the XC tables) was reached by unbounded linear
    !! extrapolation.
    subroutine test_spline2d_synthetic_node()
        use fortuno_serial, only: check => serial_check
        use spline2d, only: spline2d_t, spline2d_init, spline2d_eval, spline2d_destroy
        use lsda_constants, only: dp

        type(spline2d_t) :: spl
        integer, parameter :: nx = 4, ny_max = 3
        real(dp) :: x_grid(nx), y_grid(ny_max, nx), f_values(ny_max, nx)
        integer :: n_y_pts(nx)
        real(dp) :: f_interp
        integer :: i, j

        do i = 1, nx
            x_grid(i) = real(i, dp)
            n_y_pts(i) = ny_max
            do j = 1, ny_max
                y_grid(j, i) = real(j-1, dp) * 0.9_dp
                f_values(j, i) = x_grid(i) * x_grid(i) + y_grid(j, i)
            end do
        end do

        call spline2d_init(spl, x_grid, y_grid, f_values, n_y_pts)

        ! y below the first abscissa: the full window is used
        f_interp = spline2d_eval(spl, 0.4_dp, 0.4_dp, node0_value = 7.0_dp)
        call check(abs(f_interp - 7.0_dp) < 1.0e-12_dp, &
                   "evaluation at x = y must return the synthetic node value")

        ! y inside the grid: the window drops the rows with x < y
        f_interp = spline2d_eval(spl, 1.5_dp, 1.5_dp, node0_value = -3.0_dp)
        call check(abs(f_interp + 3.0_dp) < 1.0e-12_dp, &
                   "restricted window must still honour the synthetic node")

        ! A point between the synthetic node and the first row stays bracketed
        f_interp = spline2d_eval(spl, 0.7_dp, 0.4_dp, node0_value = 0.0_dp)
        call check(f_interp > -1.0_dp .and. f_interp < 2.0_dp, &
                   "value between the synthetic node and the first row must stay bounded")

        call spline2d_destroy(spl)

    end subroutine test_spline2d_synthetic_node

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
