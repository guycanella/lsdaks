!> Command line front end of the XC table generator.
!!
!! Solves the thermodynamic-limit Lieb-Wu integral equations over the whole
!! `(n, m)` grid for a single `U` and writes the table in the native format.
!! Refuses to write a table containing NaN or Inf and exits with status 1 in
!! that case, listing the offending grid points.
!!
!! The default output directory is the one the SCF reads, and the file name is
!! derived from `U` alone, so a run can land exactly on a table that is part of
!! the validated reference set.  An existing file is therefore never
!! overwritten unless `--force` is given.
!!
!! Usage:
!! ```
!! generate_xc_table --U <value> [--output <dir>] [--n-points N] [--m-points M]
!!                   [--m-grade G] [--n-grade-low A] [--n-grade-high B]
!!                   [--m-frac-min F] [--n-min N] [--n-k N] [--n-lambda N]
!!                   [--n-omega N] [--tol T] [--delta-n H] [--force]
!! ```
!!
!! @note The grid defaults were chosen by measuring the post-spline deviation
!!       from the reference `U = 4` table; see `bethe_tables`.  Measured wall
!!       time of the whole default grid at `U = 4`, release build on 14 cores:
!!       49 s with OpenMP enabled (`--profile release --flag -fopenmp
!!       --link-flag -fopenmp`) and 306 s (about 5 min) with a single thread.
!!       Weak coupling costs more because the `Lambda` quadrature floor rises
!!       as `U` falls.  The exact wall time depends on the adaptive `B` cutoff
!!       and the requested U; it is not a correctness criterion for a table.
program generate_xc_table_app
    use bethe_tables, only: generate_xc_table, grid_params_t, U_TABLE_MIN
    use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT
    use table_io, only: xc_table_t, write_fortran_table, count_nonfinite_entries
    use lsda_constants, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    real(dp) :: U
    character(len=256) :: output_dir, output_file, arg, val
    type(grid_params_t) :: params
    type(xc_table_t) :: table
    integer :: ierr, io_stat, nargs
    integer :: n_bad_exc, n_bad_up, n_bad_dn, i, j
    logical :: have_U, force, exists
    real(dp) :: t_start, t_end, wall
    integer(8) :: c_start, c_end, c_rate

    nargs = command_argument_count()
    output_dir = 'data/tables/fortran_native'
    params = grid_params_t()
    have_U = .false.
    force = .false.
    U = 0.0_dp

    if (nargs < 2) then
        call usage()
        stop 1
    end if

    i = 1
    do while (i <= nargs)
        call get_command_argument(i, arg)

        if (trim(arg) == '--force') then
            force = .true.
            i = i + 1
            cycle
        end if

        if (i + 1 > nargs) then
            print '(A,A)', "ERROR: missing value for option ", trim(arg)
            stop 1
        end if
        call get_command_argument(i + 1, val)

        io_stat = 0
        select case (trim(arg))
        case ('--U')
            read(val, *, iostat=io_stat) U
            have_U = .true.
        case ('--output')
            output_dir = val
        case ('--n-points')
            read(val, *, iostat=io_stat) params%n_points
        case ('--m-points')
            read(val, *, iostat=io_stat) params%m_points
        case ('--m-frac-min')
            read(val, *, iostat=io_stat) params%m_frac_min
        case ('--m-grade')
            read(val, *, iostat=io_stat) params%m_grade
        case ('--n-grade-low')
            read(val, *, iostat=io_stat) params%n_grade_low
        case ('--n-grade-high')
            read(val, *, iostat=io_stat) params%n_grade_high
        case ('--n-min')
            read(val, *, iostat=io_stat) params%n_min
        case ('--n-k')
            read(val, *, iostat=io_stat) params%quad%n_k
        case ('--n-lambda')
            read(val, *, iostat=io_stat) params%quad%n_lambda
        case ('--n-omega')
            read(val, *, iostat=io_stat) params%quad%n_omega
        case ('--tol')
            read(val, *, iostat=io_stat) params%quad%tol
        case ('--delta-n')
            read(val, *, iostat=io_stat) params%delta_n
        case default
            print '(A,A)', "ERROR: unknown option ", trim(arg)
            call usage()
            stop 1
        end select

        if (io_stat /= 0) then
            print '(A,A,A,A)', "ERROR: invalid value '", trim(val), "' for ", trim(arg)
            ierr = ERROR_INVALID_INPUT
            stop 1
        end if

        i = i + 2
    end do

    if (.not. have_U) then
        print '(A)', "ERROR: --U is required"
        call usage()
        stop 1
    end if

    if (.not. ieee_is_finite(U) .or. params%n_points < 1 .or. params%m_points < 1 .or. params%quad%n_k < 4 &
        .or. params%quad%n_lambda < 4 .or. params%quad%n_omega < 4 &
        .or. .not. ieee_is_finite(params%quad%tol) .or. params%quad%tol <= 0.0_dp &
        .or. .not. ieee_is_finite(params%delta_n) .or. params%delta_n <= 0.0_dp) then
        print '(A)', "ERROR: U must be finite; grid sizes, quadrature orders and tolerances must be positive"
        stop 1
    end if

    if (params%m_grade <= 0.0_dp .or. params%m_grade > 1.0_dp &
        .or. params%m_frac_min <= 0.0_dp .or. params%m_frac_min >= 1.0_dp &
        .or. params%n_grade_low < 1.0_dp .or. params%n_grade_high < 1.0_dp &
        .or. .not. ieee_is_finite(params%m_frac_min) .or. .not. ieee_is_finite(params%m_grade) &
        .or. .not. ieee_is_finite(params%n_grade_low) .or. .not. ieee_is_finite(params%n_grade_high) &
        .or. .not. ieee_is_finite(params%n_min) .or. .not. ieee_is_finite(params%n_max) &
        .or. params%n_min <= 0.0_dp .or. params%n_min > params%n_max &
        .or. (params%n_points > 1 .and. params%n_min >= params%n_max)) then
        print '(A)', "ERROR: --m-grade and --m-frac-min must be in (0, 1]; " // &
                     "--n-grade-low and --n-grade-high must be >= 1; " // &
                     "--n-min must satisfy 0 < n_min <= n_max (strict for multi-point grids)"
        stop 1
    end if

    ! The Lieb-Wu kernels have width U/4 in sin(k), so as U -> 0 they turn into
    ! delta functions that the fixed-order k quadrature cannot resolve: the
    ! solver itself refuses below U_QUAD_MIN = 0.5.  A *table* needs more than
    ! a point evaluation to be trustworthy - its smallest m nodes carry the
    ! m -> 0 exchange splitting - and that is only validated against the C++
    ! reference down to |U| = U_TABLE_MIN = 1.  Without this check the refusal
    ! would reach the user as a generic ERROR_INVALID_INPUT from the generator.
    !
    ! U = 0 is NOT carved out: `generate_xc_table` refuses it like any other
    ! interaction below the floor.  A U = 0 table would be identically zero and
    ! nothing consumes one - the SCF skips XC entirely there - so the honest
    ! answer is to refuse instead of writing a file of zeros.  The point
    ! evaluations `compute_E_xc` / `compute_V_xc_numerical` do still accept
    ! U = 0 and return 0; only table generation refuses.
    if (abs(U) < U_TABLE_MIN) then
        print '(A,F8.4,A)', "ERROR: |U| = ", abs(U), " is below the validated floor of the"
        print '(A,F6.2,A)', "       table generator, U_TABLE_MIN = ", U_TABLE_MIN, "."
        print '(A)', "       The spin kernels have width U/4 and become delta functions as"
        print '(A)', "       U -> 0; below U = 1 the exchange splitting V_up - V_dn of the"
        print '(A)', "       smallest m nodes is no longer accurate to a few percent and"
        print '(A)', "       there is no reference table left to validate it against."
        print '(A)', "       Use |U| >= 1.  U = 0 is refused here as well: e_xc vanishes"
        print '(A)', "       identically, so the table would be a file of zeros."
        stop 1
    end if

    ! Refuse to clobber an existing table: the default output directory is the
    ! one the SCF reads and the file name follows from U alone, so a plain
    ! `--U 4` would otherwise overwrite the C++-validated reference table.
    write(output_file, '(A,A,F0.2,A)') trim(output_dir), '/xc_table_u', U, '.dat'
    inquire(file=trim(output_file), exist=exists)
    if (exists .and. .not. force) then
        print '(A,A)', "ERROR: refusing to overwrite existing file ", trim(output_file)
        print '(A)', "       Files under data/tables/fortran_native are the reference tables"
        print '(A)', "       validated against the C++ implementation."
        print '(A)', "       Write elsewhere with --output <dir>, or pass --force if you"
        print '(A)', "       really mean to replace the reference."
        stop 1
    end if

    print '(A)', "=========================================="
    print '(A)', "  XC Table Generator"
    print '(A)', "=========================================="
    print '(A,F0.2)', "  U value:     ", U
    print '(A,A)', "  Output dir:  ", trim(output_dir)
    print '(A)', ""
    print '(A)', "Grid parameters:"
    print '(A,F0.4,A,F0.4)', "  n range:  ", params%n_min, " to ", params%n_max
    print '(A,I0,A,I0,A,I0,A)', "  Grid:     ", params%n_points, " x ", &
                                params%m_points, " = ", &
                                params%n_points * params%m_points, " points"
    print '(A,F6.3,A,F6.3,A)', "  n grading: low end ", params%n_grade_low, &
                               ", high end ", params%n_grade_high, " (1 = uniform)"
    print '(A,F6.3,A,ES9.2)', "  m grading: ", params%m_grade, &
                              " (1 = uniform), smallest m/n = ", params%m_frac_min
    print '(A)', ""
    print '(A)', "Numerics (thermodynamic-limit Lieb-Wu integral equations):"
    print '(A,I0,A,I0,A,I0)', "  Quadrature: n_k = ", params%quad%n_k, &
        ", n_lambda = ", params%quad%n_lambda, ", n_omega = ", params%quad%n_omega
    print '(A,ES9.2)', "  Inversion tolerance: ", params%quad%tol
    print '(A,ES9.2)', "  V_xc finite-difference step: ", params%delta_n
    print '(A)', ""

    print '(A)', "Starting table generation..."
    print '(A)', ""

    call cpu_time(t_start)
    call system_clock(c_start, c_rate)
    call generate_xc_table(U, params, table, ierr)
    call system_clock(c_end)
    call cpu_time(t_end)

    if (ierr /= ERROR_SUCCESS) then
        print *, "ERROR: Table generation failed!"
        stop 1
    end if

    ! Wall time, not CPU time: the grid loop is an OpenMP parallel do, so
    ! cpu_time sums over all threads and says nothing about how long the user
    ! actually waited.
    wall = real(c_end - c_start, dp) / real(max(c_rate, 1_8), dp)
    print '(A,F8.2,A,F8.2,A)', "Elapsed: ", wall, " s wall, ", t_end - t_start, " s CPU"

    ! write_fortran_table also refuses non-finite tables; the check is done
    ! here first so the failing grid points can be listed for the user.
    call count_nonfinite_entries(table, n_bad_exc, n_bad_up, n_bad_dn)
    if (n_bad_exc + n_bad_up + n_bad_dn > 0) then
        print '(A)', ""
        print '(A)', "ERROR: Generated table contains non-finite entries (NaN or Inf); nothing was written."
        print '(A,I0,A,I0,A,I0)', "       non-finite count: exc = ", n_bad_exc, &
            ", Vxc_up = ", n_bad_up, ", Vxc_down = ", n_bad_dn
        print '(A)', "       Failing grid points (n, m):"
        do i = 1, table%n_points_n
            do j = 1, table%n_points_m
                if (.not. (ieee_is_finite(table%exc(j, i)) .and. ieee_is_finite(table%vxc_up(j, i)) &
                    .and. ieee_is_finite(table%vxc_down(j, i)))) then
                    print '(A,F8.4,A,F8.4)', "         n = ", table%n_grid(i), &
                        "   m = ", table%m_grid(j, i)
                end if
            end do
        end do
        stop 1
    end if

    print '(A)', ""
    print '(A)', "Saving table to disk..."

    call write_fortran_table(output_file, table, ierr)

    if (ierr /= ERROR_SUCCESS) then
        print *, "ERROR: Failed to write table file!"
        stop 1
    end if

    print '(A)', ""
    print '(A)', "=========================================="
    print '(A)', "Table written (all entries finite)."
    print '(A)', "=========================================="
    print '(A,A)', "Table saved to: ", trim(output_file)
    print '(A)', ""

contains

    subroutine usage()
        print '(A)', "Usage: generate_xc_table --U <value> [options]"
        print '(A)', ""
        print '(A)', "Options:"
        print '(A)', "  --U <value>        Hubbard interaction (required, may be negative)"
        print '(A)', "  --output <dir>     Output directory"
        print '(A)', "  --n-points <int>   Density grid points"
        print '(A)', "  --m-points <int>   Magnetization grid points per density"
        print '(A)', "  --m-grade <real>   m axis grading exponent in (0, 1]; 1 = uniform"
        print '(A)', "  --m-frac-min <real>  Smallest non-zero m/n of the m axis"
        print '(A)', "  --n-grade-low <real>   n axis grading at n_min; >= 1, 1 = uniform"
        print '(A)', "  --n-grade-high <real>  n axis grading at n_max; >= 1, 1 = uniform"
        print '(A)', "  --n-min <real>     Lowest density of the grid"
        print '(A)', "  --n-k <int>        Gauss-Legendre nodes on [0, Q]"
        print '(A)', "  --n-lambda <int>   Gauss-Legendre nodes per Lambda panel"
        print '(A)', "  --n-omega <int>    Gauss-Legendre nodes per omega panel"
        print '(A)', "  --tol <real>       Density residual tolerance of the (Q, B) inversion"
        print '(A)', "  --delta-n <real>   Finite-difference step of V_xc"
        print '(A)', "  --force            Overwrite an existing table file"
        print '(A)', ""
        print '(A)', "Example: fpm run generate_xc_table --profile release --flag -fopenmp" // &
                     " -- --U 5.5"
    end subroutine usage

end program generate_xc_table_app
