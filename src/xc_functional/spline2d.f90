!> Two-dimensional cubic spline interpolation on a row-wise irregular grid
!!
!! The grid is a stack of rows: row `i` sits at abscissa `x(i)` and carries its
!! own y-grid `y(1:n_y(i), i)`. This is the layout of the Bethe Ansatz XC
!! tables, where the magnetization grid of each density row spans m ∈ [0, n]
!! and therefore changes from row to row.
!!
!! Interpolation is **cubic in both directions** (it used to be cubic in y and
!! merely linear in x):
!!  1. `spline2d_init` pre-computes, for every row, the second derivatives of
!!     the cubic spline in the y direction.
!!  2. `spline2d_eval` evaluates those row splines at the requested y, which
!!     produces a 1D data set (x(i), f_i), and then builds a *clamped* cubic
!!     spline through it and evaluates it at the requested x.
!!
!! The x-direction spline is assembled at every call, exactly as in the C++
!! reference (`original/spline2D.cc`, functions `exc_value`, `Vxc_up_value`
!! and `Vxc_dn_value`). That reference also restricts the x window to the rows
!! that actually bracket the requested y and prepends one synthetic node; both
!! features are available through the optional arguments of `spline2d_eval`.
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
    private :: eval_row

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
    !! The boundary condition of the row splines is selectable. The default,
    !! `'natural'`, sets the second derivative to zero at both ends of every
    !! row. The C++ reference (`original/spline2D.cc`, `build_spline_data`)
    !! instead clamps them: the two V_xc tables get the secants of the first and
    !! last interval, and the e_xc table gets the physically known derivatives
    !! (zero at m = 0 by spin symmetry, and the analytic slope of the fully
    !! polarized line at m = n). Passing `row_bc = 'clamped'` reproduces that;
    !! `dy_first` / `dy_last` override the secant of the corresponding end.
    !!
    !! @param[out] spl      Spline object to initialize
    !! @param[in]  x_grid   Grid points in x direction (n)
    !! @param[in]  y_grid   Grid points in y direction (m), shape (n_y_max, n_x)
    !! @param[in]  f_values Function values on grid, shape (n_y_max, n_x)
    !! @param[in]  n_y_pts  Number of valid y points for each x, shape (n_x)
    !! @param[in]  row_bc   Optional row boundary condition, `'natural'`
    !!                      (default) or `'clamped'`
    !! @param[in]  dy_first Optional prescribed df/dy at the FIRST y point of
    !!                      each row, shape (n_x); only used with
    !!                      `row_bc = 'clamped'` (default: first secant)
    !! @param[in]  dy_last  Optional prescribed df/dy at the LAST y point of
    !!                      each row, shape (n_x); only used with
    !!                      `row_bc = 'clamped'` (default: last secant)
    subroutine spline2d_init(spl, x_grid, y_grid, f_values, n_y_pts, row_bc, dy_first, dy_last)
        type(spline2d_t), intent(out) :: spl
        real(dp), intent(in) :: x_grid(:), y_grid(:,:), f_values(:,:)
        integer, intent(in) :: n_y_pts(:)
        character(len=*), intent(in), optional :: row_bc
        real(dp), intent(in), optional :: dy_first(:), dy_last(:)

        !> Two y nodes closer than this are treated as coincident
        real(dp), parameter :: Y_DEGENERATE_TOL = 1.0e-15_dp

        integer :: i, nx, ny_max, ny
        character(len=32) :: bc
        real(dp) :: dy0, dyn

        bc = 'natural'
        if (present(row_bc)) bc = row_bc

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

                    if (trim(bc) == 'clamped') then
                        ! Default: the secants of the extreme intervals; a
                        ! degenerate interval (padded row) falls back to a flat
                        ! end instead of dividing by zero.
                        dy0 = 0.0_dp
                        if (abs(y_temp(1) - y_temp(0)) > Y_DEGENERATE_TOL) &
                            dy0 = (f_temp(1) - f_temp(0)) / (y_temp(1) - y_temp(0))

                        dyn = 0.0_dp
                        if (abs(y_temp(ny-1) - y_temp(ny-2)) > Y_DEGENERATE_TOL) &
                            dyn = (f_temp(ny-1) - f_temp(ny-2)) / (y_temp(ny-1) - y_temp(ny-2))

                        if (present(dy_first)) dy0 = dy_first(i)
                        if (present(dy_last)) dyn = dy_last(i)

                        call spline1d_coeff(y_temp, f_temp, ny-1, d2y_temp, 'clamped', dy0, dyn)
                    else
                        call spline1d_coeff(y_temp, f_temp, ny-1, d2y_temp, 'natural', 0.0_dp, 0.0_dp)
                    end if

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

    !> Evaluate the y-direction spline of one row
    !!
    !! Evaluates the pre-computed cubic spline of row `i_x` at the point `y`.
    !! Outside the row's y range the polynomial of the extreme interval is
    !! continued (same behaviour as `spline3` in the C++ reference).
    !!
    !! @param[in] spl  Initialized spline object
    !! @param[in] i_x  Row index (1 ≤ i_x ≤ spl%n_x)
    !! @param[in] y    Y-coordinate
    !! @return         Interpolated value f(x(i_x), y)
    function eval_row(spl, i_x, y) result(f_row)
        type(spline2d_t), intent(in) :: spl
        integer, intent(in) :: i_x
        real(dp), intent(in) :: y
        real(dp) :: f_row

        integer :: ny

        ny = spl%n_y(i_x)
        f_row = spline1d_eval(spl%y(1:ny, i_x), spl%f(1:ny, i_x), &
                              spl%d2f_dy2(1:ny, i_x), ny - 1, y)
    end function eval_row

    !> Evaluate 2D spline at point (x, y) with cubic interpolation in BOTH directions
    !!
    !! Algorithm:
    !! 1. Select the window of rows used in the x direction (see below).
    !! 2. Evaluate the y-direction spline of every row of the window at `y`,
    !!    producing the 1D data set (x_i, f_i).
    !! 3. Build a *clamped* cubic spline through that data set and evaluate it
    !!    at `x`.
    !!
    !! Without the optional arguments the window is the whole set of rows and
    !! both prescribed end derivatives are the corresponding secants; this is
    !! the general-purpose bicubic interpolation.
    !!
    !! With `node0_value` present the routine reproduces literally the scheme
    !! of the C++ reference for the XC tables:
    !! - the window starts at the last row whose abscissa is still ≥ `y`
    !!   (found by the same linear scan as the original), and runs to the last
    !!   tabulated row; the window therefore *grows* as `y` decreases;
    !! - a synthetic node (x = y, f = `node0_value`) is prepended. In the XC
    !!   tables x is the total density n and y is the magnetization m, so
    !!   x = y is the fully polarized line n_dn = 0, whose value is known
    !!   analytically. Because n ≥ |m| always holds, this node also removes
    !!   any extrapolation below the first tabulated density;
    !! - if the window contains a single tabulated row the original falls back
    !!   to linear interpolation; pass `allow_linear_branch = .false.` to
    !!   suppress that fallback (the C++ `Vxc_dn_value` has no such branch).
    !!
    !! @param[in] spl                 Initialized spline object
    !! @param[in] x                   X-coordinate (density n)
    !! @param[in] y                   Y-coordinate (magnetization m)
    !! @param[in] node0_value         Optional value of the synthetic node placed
    !!                                at x = y; its presence also activates the
    !!                                restricted row window
    !! @param[in] dfdx_first          Optional prescribed df/dx at the first node
    !!                                (default: the secant of the first interval)
    !! @param[in] allow_linear_branch Optional flag (default .true.) enabling the
    !!                                linear fallback for a single-row window
    !! @return                        Interpolated value f(x, y)
    function spline2d_eval(spl, x, y, node0_value, dfdx_first, allow_linear_branch) result(f_interp)
        type(spline2d_t), intent(in) :: spl
        real(dp), intent(in) :: x, y
        real(dp), intent(in), optional :: node0_value
        real(dp), intent(in), optional :: dfdx_first
        logical, intent(in), optional :: allow_linear_branch
        real(dp) :: f_interp

        !> Two abscissae closer than this are treated as coincident
        real(dp), parameter :: X_DEGENERATE_TOL = 1.0e-15_dp

        real(dp), allocatable :: x_loc(:), f_loc(:), d2_loc(:)
        real(dp) :: dy_ini, dy_fim
        integer :: i, i_in, i_first, num_n, offset
        logical :: linear_ok, with_node0

        if (.not. spl%initialized) then
            print *, "ERROR: spline2d not initialized!"
            f_interp = 0.0_dp
            return
        end if

        linear_ok = .true.
        if (present(allow_linear_branch)) linear_ok = allow_linear_branch

        with_node0 = present(node0_value)

        if (with_node0) then
            ! Same scan as the C++ reference: advance to the first row with
            ! x(i) >= y, then step one row back.
            i_in = 1
            do while (spl%x(i_in) < y .and. i_in < spl%n_x)
                i_in = i_in + 1
            end do
            i_in = i_in - 1
            i_first = i_in + 1

            ! If `y` coincides with a tabulated abscissa the synthetic node
            ! would duplicate the first row and every secant would divide by
            ! zero (the C++ reference divides by zero here). Drop the synthetic
            ! node instead: it carries no information in that limit.
            if (abs(spl%x(i_first) - y) < X_DEGENERATE_TOL) with_node0 = .false.
        else
            i_first = 1
        end if

        if (with_node0) then
            num_n = spl%n_x - i_first + 1
            offset = 1
        else
            num_n = spl%n_x - i_first
            offset = 0
        end if

        allocate(x_loc(0:num_n), f_loc(0:num_n), d2_loc(0:num_n))

        if (with_node0) then
            x_loc(0) = y
            f_loc(0) = node0_value
        end if

        do i = offset, num_n
            x_loc(i) = spl%x(i_first + i - offset)
            f_loc(i) = eval_row(spl, i_first + i - offset, y)
        end do

        if (num_n < 1) then
            ! Single node: nothing to interpolate in x
            f_interp = f_loc(0)
        else if (num_n == 1 .and. linear_ok) then
            f_interp = f_loc(0) + (f_loc(1) - f_loc(0)) / (x_loc(1) - x_loc(0)) * (x - x_loc(0))
        else
            dy_ini = (f_loc(1) - f_loc(0)) / (x_loc(1) - x_loc(0))
            if (present(dfdx_first)) dy_ini = dfdx_first
            dy_fim = (f_loc(num_n) - f_loc(num_n-1)) / (x_loc(num_n) - x_loc(num_n-1))

            call spline1d_coeff(x_loc, f_loc, num_n, d2_loc, 'clamped', dy_ini, dy_fim)
            f_interp = spline1d_eval(x_loc, f_loc, d2_loc, num_n, x)
        end if

        deallocate(x_loc, f_loc, d2_loc)
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
