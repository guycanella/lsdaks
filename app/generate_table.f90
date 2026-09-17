program generate_xc_table_app
    use bethe_tables, only: generate_xc_table, grid_params_t
    use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT
    use table_io, only: xc_table_t, write_fortran_table, count_nonfinite_entries
    use lsda_constants, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    real(dp) :: U
    character(len=256) :: output_dir, output_file, arg
    type(grid_params_t) :: params
    type(xc_table_t) :: table
    integer :: ierr, io_stat, nargs
    integer :: n_bad_exc, n_bad_up, n_bad_dn, i, j
    
    nargs = command_argument_count()
    
    if (nargs < 2) then
        print *, "Usage: generate_xc_table --U <value> [--output <dir>]"
        print *, "Example: fpm run generate_xc_table -- --U 5.5"
        stop 1
    end if
    
    call get_command_argument(1, arg)
    if (trim(arg) /= '--U') then
        print *, "ERROR: First argument must be --U"
        stop 1
    end if
    
    call get_command_argument(2, arg)
    read(arg, *, iostat=io_stat) U
    if (io_stat /= 0) then
        print *, "ERROR: Invalid U value:", trim(arg)
        ierr = ERROR_INVALID_INPUT
        stop 1
    end if
    
    output_dir = 'data/tables/fortran_native'
    if (nargs >= 4) then
        call get_command_argument(3, arg)
        if (trim(arg) == '--output') then
            call get_command_argument(4, output_dir)
        end if
    end if
    
    print '(A)', "=========================================="
    print '(A)', "  XC Table Generator"
    print '(A)', "=========================================="
    print '(A,F0.2)', "  U value:     ", U
    print '(A,A)', "  Output dir:  ", trim(output_dir)
    print '(A)', ""
    
    params = grid_params_t()
    
    print '(A)', "Grid parameters:"
    print '(A,F0.2,A,F0.2)', "  n range:  ", params%n_min, " to ", params%n_max
    print '(A,I0,A,I0,A,I0)', "  Grid:     ", params%n_points, " x ", &
                               params%m_points, " = ", &
                               params%n_points * params%m_points, " points"
    print '(A,I0)', "  Sys size: L = ", params%L
    print '(A)', ""
    print '(A)', "⚠️  EXPERIMENTAL: finite-L Bethe Ansatz generator."
    print '(A)', "    The Newton solver still fails to converge on part of the grid"
    print '(A)', "    and the finite-L V_xc carries O(1/L) error. Tables produced"
    print '(A)', "    here must NOT be used in production SCF runs; use the"
    print '(A)', "    reference tables in data/tables/fortran_native instead."
    print '(A)', "    See NEXT_STEPS_REPORT.md, phase 4.5 (T21)."
    print '(A)', ""

    print '(A)', "Starting table generation..."
    print '(A)', ""

    call generate_xc_table(U, params, table, ierr)

    if (ierr /= ERROR_SUCCESS) then
        print *, "ERROR: Table generation failed!"
        stop 1
    end if

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

    ! Construct full filename
    write(output_file, '(A,A,F0.2,A)') trim(output_dir), '/xc_table_u', U, '.dat'

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
    print '(A)', "Reminder: this generator is EXPERIMENTAL (finite-L, see phase 4.5)."
    print '(A)', "Validate against a reference table before using it in an SCF run."
    print '(A)', ""
    
end program generate_xc_table_app
