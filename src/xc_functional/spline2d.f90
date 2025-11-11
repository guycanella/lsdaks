module spline2d
    use lsda_constants, only: dp
    implicit none
    private

    type, public :: spline2d_t
          integer :: n_x                        !< Number of points in x (n) direction
          integer :: n_y_max                    !< Maximum number of points in y (m) direction
          real(dp), allocatable :: x(:)         !< Grid points in x direction: x(n_x)
          real(dp), allocatable :: y(:,:)       !< Grid points in y direction: y(n_y, n_x) - varies with x!
          real(dp), allocatable :: f(:,:)       !< Function values: f(n_y, n_x)
          integer, allocatable :: n_y(:)        !< Number of y points for each x: n_y(n_x)

          real(dp), allocatable :: d2f_dy2(:,:) !< Second derivatives in y: d2f_dy2(n_y, n_x)

          logical :: initialized = .false.      !< Initialization flag
    end type spline2d_t

    public :: spline2d_init
    public :: spline2d_eval
    public :: spline2d_destroy

    private :: spline1d_coeff
    private :: spline1d_eval
    private :: find_interval

contains
    !> Compute cubic spline coefficients (second derivatives) for 1D data
    !!
    !! Solves tridiagonal system for cubic spline with specified boundary conditions.
    !!
    !! Boundary conditions:
    !! - 'natural': second derivative = 0 at endpoints (natural spline)
    !! - 'clamped': first derivative specified at endpoints (clamped spline)
    !!
    !! Based on the C++ code "gera_dy2_alternativa" from the original implementation.
    !!
    !! @param[in]  x       Grid points, monotonically increasing, size n+1 (0 to n)
    !! @param[in]  y       Function values at grid points, size n+1
    !! @param[in]  n       Number of intervals (n+1 points total)
    !! @param[out] d2y     Second derivatives at grid points, size n+1
    !! @param[in]  bc_type Boundary condition type: 'natural' or 'clamped'
    !! @param[in]  dy0     First derivative at x(0) (used if bc_type='clamped')
    !! @param[in]  dyn     First derivative at x(n) (used if bc_type='clamped')
    subroutine spline1d_coeff(x, y, n, d2y, bc_type, dy0, dyn)
        real(dp), intent(in) :: x(0:), y(0:), dy0, dyn
        integer, intent(in) :: n
        real(dp), intent(out) :: d2y(0:)
        character(len=*), intent(in) :: bc_type

        real(dp), allocatable :: a(:), b(:)
        real(dp) :: h0, h1
        integer :: i

        allocate(a(0:n), b(0:n))

        if (trim(bc_type) == 'clamped') then
            a(0) = (x(1) - x(0)) / 3.0_dp
            b(0) = (y(1) - y(0)) / (x(1) - x(0)) - dy0

            do i = 1, n-1
                a(i) = (x(i+1) - x(i-1)) / 3.0_dp
                b(i) = (y(i+1) - y(i)) / (x(i+1) - x(i)) - &
                        (y(i) - y(i-1)) / (x(i) - x(i-1))
            end do

            a(n) = (x(n) - x(n-1)) / 3.0_dp
            b(n) = dyn - (y(n) - y(n-1)) / (x(n) - x(n-1))

        else
            a(0) = 1.0_dp
            b(0) = 0.0_dp

            do i = 1, n-1
                a(i) = (x(i+1) - x(i-1)) / 3.0_dp
                b(i) = (y(i+1) - y(i)) / (x(i+1) - x(i)) - &
                        (y(i) - y(i-1)) / (x(i) - x(i-1))
            end do

            a(n) = 1.0_dp
            b(n) = 0.0_dp
        end if

        ! Forward elimination (Thomas algorithm for tridiagonal system)
        do i = 1, n-1
            h0 = x(i) - x(i-1)
            h1 = x(i+1) - x(i)

            if (trim(bc_type) == 'clamped') then
                a(i) = a(i) - h0 * h0 / a(i-1) / 36.0_dp
                b(i) = b(i) - b(i-1) / a(i-1) * h0 / 6.0_dp
            else
                if (i == 1) then
                    cycle
                end if

                a(i) = a(i) - h0 * h0 / a(i-1) / 36.0_dp
                b(i) = b(i) - b(i-1) / a(i-1) * h0 / 6.0_dp
            end if
        end do

        if (trim(bc_type) == 'clamped') then
            h0 = x(n) - x(n-1)
            a(n) = h0 / 3.0_dp * (1.0_dp - h0 / a(n-1) / 12.0_dp)
            b(n) = b(n) - b(n-1) / a(n-1) * h0 / 6.0_dp
        end if

        ! Back substitution
        d2y(n) = b(n) / a(n)

        do i = n-1, 1, -1
            h1 = (x(i+1) - x(i)) / 6.0_dp
            d2y(i) = (b(i) - h1 * d2y(i+1)) / a(i)
        end do

        if (trim(bc_type) == 'clamped') then
            d2y(0) = b(0) / a(0) - d2y(1) / 2.0_dp
        else
            d2y(0) = 0.0_dp
        end if

        deallocate(a, b)
    end subroutine spline1d_coeff

    !> Initialize 2D spline from grid data
    !!
    !! Constructs cubic spline interpolation for irregular 2D grid.
    !! For each fixed x_i, computes spline coefficients in y direction.
    !!
    !! @param[out] spl      Spline object to initialize
    !! @param[in]  x_grid   Grid points in x direction (n)
    !! @param[in]  y_grid   Grid points in y direction (m), shape (n_y_max, n_x)
    !! @param[in]  f_values Function values on grid, shape (n_y_max, n_x)
    !! @param[in]  n_y_pts  Number of valid y points for each x, shape (n_x)
    subroutine spline2d_init(spl, x_grid, y_grid, f_values, n_y_pts)
        type(spline2d_t), intent(out) :: spl
        real(dp), intent(in) :: x_grid(:), y_grid(:,:), f_values(:,:)
        integer, intent(in) :: n_y_pts(:)

        integer :: i, nx, ny_max, ny

        nx = size(x_grid)
        ny_max = size(y_grid, 1)

        spl%n_x = nx
        spl%n_y_max = ny_max

        allocate(spl%x(nx))
        allocate(spl%y(ny_max, nx))
        allocate(spl%f(ny_max, nx))
        allocate(spl%n_y(nx))
        allocate(spl%d2f_dy2(ny_max, nx))

        spl%x = x_grid
        spl%y = y_grid
        spl%f = f_values
        spl%n_y = n_y_pts

        do i = 1, nx
            ny = spl%n_y(i)
            if (ny > 1) then
                ! Create temporary 0-indexed arrays for spline1d_coeff
                block
                    real(dp) :: y_temp(0:ny-1), f_temp(0:ny-1), d2y_temp(0:ny-1)

                    y_temp = spl%y(1:ny, i)
                    f_temp = spl%f(1:ny, i)

                    call spline1d_coeff(y_temp, f_temp, ny-1, d2y_temp, 'natural', 0.0_dp, 0.0_dp)

                    spl%d2f_dy2(1:ny, i) = d2y_temp
                end block
            else
                spl%d2f_dy2(1, i) = 0.0_dp
            end if
        end do

        spl%initialized = .true.
    end subroutine spline2d_init

    !> Find interval containing x in monotonic array
    !!
    !! Binary search to find i such that x_grid(i) <= x < x_grid(i+1)
    !! Array is assumed to be 0-indexed: x_grid(0), x_grid(1), ..., x_grid(n-1)
    !!
    !! @param[in] x_grid Monotonically increasing array (0 to n-1)
    !! @param[in] n      Array size
    !! @param[in] x      Point to locate
    !! @return           Index i of left endpoint (0 to n-2)
    function find_interval(x_grid, n, x) result(i)
        real(dp), intent(in) :: x_grid(0:), x
        integer, intent(in) :: n
        integer :: left, right, mid, i

        if (x <= x_grid(0)) then
            i = 0
            return
        end if

        if (x >= x_grid(n-1)) then
            i = n - 2
            return
        end if

        left = 0
        right = n - 1

        do while (right - left > 1)
            mid = (left + right) / 2
            if (x < x_grid(mid)) then
                right = mid
            else
                left = mid
            end if
        end do

        i = left
    end function find_interval

    !> Evaluate cubic spline at point x given coefficients
    !!
    !! Uses the cubic spline formula with second derivatives.
    !!
    !! @param[in] x_grid Grid points (0 to n)
    !! @param[in] y_grid Function values (0 to n)
    !! @param[in] d2y    Second derivatives (from spline1d_coeff)
    !! @param[in] n      Number of intervals (n+1 points)
    !! @param[in] x      Point to evaluate
    !! @return           Interpolated value
    function spline1d_eval(x_grid, y_grid, d2y, n, x) result(y_interp)
        real(dp), intent(in) :: x_grid(0:), y_grid(0:), d2y(0:), x
        integer, intent(in) :: n
        real(dp) :: y_interp, h, a, b

        integer :: i

        if (n == 0) then
            y_interp = y_grid(0)
            return
        end if

        i = find_interval(x_grid, n+1, x)

        if (i < 0) i = 0
        if (i >= n) i = n - 1

        h = x_grid(i+1) - x_grid(i)

        if (abs(h) < 1.0e-15_dp) then
            y_interp = y_grid(i)
            return
        end if

        a = (x_grid(i+1) - x) / h
        b = 1.0_dp - a

        y_interp = a * y_grid(i) + b * y_grid(i+1) + &
                    ((a**3 - a) * h**2 / 6.0_dp) * d2y(i) + &
                    ((b**3 - b) * h**2 / 6.0_dp) * d2y(i+1)
    end function spline1d_eval

    !> Evaluate 2D spline at point (x, y)
    !!
    !! Uses separable interpolation:
    !! 1. Find x interval: x_i <= x <= x_{i+1}
    !! 2. Evaluate 1D spline at (x_i, y) and (x_{i+1}, y)
    !! 3. Linearly interpolate between these values in x direction
    !!
    !! @param[in] spl Initialized spline object
    !! @param[in] x   X-coordinate (density n)
    !! @param[in] y   Y-coordinate (magnetization m)
    !! @return        Interpolated value f(x, y)
    function spline2d_eval(spl, x, y) result(f_interp)
        type(spline2d_t), intent(in) :: spl
        real(dp), intent(in) :: x, y
        real(dp) :: f_interp, f_left, f_right, t

        integer :: i_x

        if (.not. spl%initialized) then
            print *, "ERROR: spline2d not initialized!"
            f_interp = 0.0_dp
            return
        end if

        i_x = find_interval(spl%x, spl%n_x, x)

        ! Convert 0-based index to 1-based for Fortran arrays
        i_x = i_x + 1

        if (spl%n_x == 1) then
            ! Single x point - evaluate 1D spline in y direction
            block
                integer :: ny1
                real(dp) :: y_tmp(0:spl%n_y(1)-1), f_tmp(0:spl%n_y(1)-1), d2y_tmp(0:spl%n_y(1)-1)
                ny1 = spl%n_y(1)
                y_tmp = spl%y(1:ny1, 1)
                f_tmp = spl%f(1:ny1, 1)
                d2y_tmp = spl%d2f_dy2(1:ny1, 1)
                f_interp = spline1d_eval(y_tmp, f_tmp, d2y_tmp, ny1-1, y)
            end block
            return
        end if

        ! Clamp i_x to valid range
        if (i_x < 1) i_x = 1
        if (i_x >= spl%n_x) i_x = spl%n_x - 1

        ! Evaluate at left endpoint (i_x)
        block
            integer :: ny_left
            real(dp), allocatable :: y_tmp(:), f_tmp(:), d2y_tmp(:)
            ny_left = spl%n_y(i_x)
            allocate(y_tmp(0:ny_left-1), f_tmp(0:ny_left-1), d2y_tmp(0:ny_left-1))
            y_tmp = spl%y(1:ny_left, i_x)
            f_tmp = spl%f(1:ny_left, i_x)
            d2y_tmp = spl%d2f_dy2(1:ny_left, i_x)
            f_left = spline1d_eval(y_tmp, f_tmp, d2y_tmp, ny_left-1, y)
            deallocate(y_tmp, f_tmp, d2y_tmp)
        end block

        ! Evaluate at right endpoint (i_x+1)
        block
            integer :: ny_right
            real(dp), allocatable :: y_tmp(:), f_tmp(:), d2y_tmp(:)
            ny_right = spl%n_y(i_x+1)
            allocate(y_tmp(0:ny_right-1), f_tmp(0:ny_right-1), d2y_tmp(0:ny_right-1))
            y_tmp = spl%y(1:ny_right, i_x+1)
            f_tmp = spl%f(1:ny_right, i_x+1)
            d2y_tmp = spl%d2f_dy2(1:ny_right, i_x+1)
            f_right = spline1d_eval(y_tmp, f_tmp, d2y_tmp, ny_right-1, y)
            deallocate(y_tmp, f_tmp, d2y_tmp)
        end block

        t = (x - spl%x(i_x)) / (spl%x(i_x+1) - spl%x(i_x))
        f_interp = (1.0_dp - t) * f_left + t * f_right
    end function spline2d_eval

    !> Clean up spline object
    subroutine spline2d_destroy(spl)
        type(spline2d_t), intent(inout) :: spl

        if (allocated(spl%x)) deallocate(spl%x)
        if (allocated(spl%y)) deallocate(spl%y)
        if (allocated(spl%f)) deallocate(spl%f)
        if (allocated(spl%n_y)) deallocate(spl%n_y)
        if (allocated(spl%d2f_dy2)) deallocate(spl%d2f_dy2)

        spl%initialized = .false.
    end subroutine spline2d_destroy
end module spline2d