program generate_xc_table_app
    use bethe_tables, only: generate_xc_table, grid_params_t
    use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT
    use table_io, only: xc_table_t, write_fortran_table
    use lsda_constants, only: dp
    implicit none

    real(dp) :: U
    character(len=256) :: output_dir, output_file, arg
    type(grid_params_t) :: params
    type(xc_table_t) :: table
    integer :: ierr, io_stat, nargs
    
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
    print '(A)', "  XC Table Generator (OpenMP Parallel)"
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
    print '(A)', "⚠️  WARNING: This operation is computationally expensive!"
    print '(A)', "    Estimated time: 10-60 minutes (depends on CPU cores)"
    print '(A)', "    OpenMP parallelization is ENABLED"
    print '(A)', ""

    print '(A)', "Starting table generation..."
    print '(A)', ""
    
    call generate_xc_table(U, params, table, ierr)
    
    if (ierr /= ERROR_SUCCESS) then
        print *, "ERROR: Table generation failed!"
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
    print '(A)', "✅ SUCCESS!"
    print '(A)', "=========================================="
    print '(A,A)', "Table saved to: ", trim(output_file)
    print '(A)', ""
    print '(A)', "You can now run simulations with:"
    print '(A,F0.2)', "  fpm run lsdaks -- --U ", U
    print '(A)', ""
    
end program generate_xc_table_app