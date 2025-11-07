!> Utility program to convert C++ ASCII tables to Fortran binary format
!!
!! This program reads all legacy C++ XC tables from a directory and converts
!! them to the native Fortran binary format for faster I/O.
!!
!! Usage:
!!   fpm run convert_tables -- <input_dir> <output_dir>
!!
!! Example:
!!   fpm run convert_tables -- data/tables/cpp_legacy data/tables/fortran_native

program convert_tables
    use iso_fortran_env, only: output_unit, error_unit
    use lsda_constants, only: dp
    use table_io
    implicit none

    character(len=256) :: input_dir, output_dir
    character(len=256) :: input_file, output_file
    type(xc_table_t) :: table
    integer :: status, n_converted, n_failed
    integer :: i

    ! List of all available U values
    real(dp), parameter :: U_VALUES(*) = [ &
        0.90_dp, 1.00_dp, 1.10_dp, 2.00_dp, 3.00_dp, 4.00_dp, 4.10_dp, &
        5.00_dp, 5.90_dp, 6.00_dp, 6.10_dp, 6.90_dp, 7.00_dp, 7.10_dp, &
        7.90_dp, 8.00_dp, 8.10_dp, 8.90_dp, 9.00_dp, 9.10_dp, 10.00_dp, &
        12.00_dp, 14.00_dp, 16.00_dp, 18.00_dp, 20.00_dp ]

    integer, parameter :: N_TABLES = size(U_VALUES)

    ! Parse command line arguments
    call parse_arguments(input_dir, output_dir, status)
    if (status /= 0) then
        call print_usage()
        stop 1
    end if

    ! Print header
    write(output_unit, '(A)') ""
    write(output_unit, '(A)') "========================================="
    write(output_unit, '(A)') "  XC Table Conversion Utility"
    write(output_unit, '(A)') "  C++ ASCII → Fortran Binary"
    write(output_unit, '(A)') "========================================="
    write(output_unit, '(A,A)') "  Input directory:  ", trim(input_dir)
    write(output_unit, '(A,A)') "  Output directory: ", trim(output_dir)
    write(output_unit, '(A,I0)') "  Number of tables: ", N_TABLES
    write(output_unit, '(A)') "========================================="
    write(output_unit, '(A)') ""

    ! Convert each table
    n_converted = 0
    n_failed = 0

    do i = 1, N_TABLES
        ! Construct filenames
        call construct_cpp_filename(input_dir, U_VALUES(i), input_file)
        call construct_fortran_filename(output_dir, U_VALUES(i), output_file)

        write(output_unit, '(A,F6.2,A)', advance='no') "  Converting U = ", U_VALUES(i), " ... "

        ! Read C++ table
        call read_cpp_table(input_file, table, status)
        if (status /= 0) then
            write(output_unit, '(A)') "FAILED (read error)"
            write(error_unit, '(A,A)') "    ERROR: Could not read file: ", trim(input_file)
            n_failed = n_failed + 1
            cycle
        end if

        ! Write Fortran binary table
        call write_fortran_table(output_file, table, status)
        if (status /= 0) then
            write(output_unit, '(A)') "FAILED (write error)"
            write(error_unit, '(A,A)') "    ERROR: Could not write file: ", trim(output_file)
            n_failed = n_failed + 1
            call deallocate_table(table)
            cycle
        end if

        ! Success
        write(output_unit, '(A)') "OK"
        n_converted = n_converted + 1

        ! Clean up
        call deallocate_table(table)
    end do

    ! Print summary
    write(output_unit, '(A)') ""
    write(output_unit, '(A)') "========================================="
    write(output_unit, '(A)') "  Conversion Summary"
    write(output_unit, '(A)') "========================================="
    write(output_unit, '(A,I0)') "  Successfully converted: ", n_converted
    write(output_unit, '(A,I0)') "  Failed:                 ", n_failed
    write(output_unit, '(A)') "========================================="
    write(output_unit, '(A)') ""

    if (n_failed > 0) then
        write(error_unit, '(A)') "WARNING: Some conversions failed. Check error messages above."
        stop 1
    else
        write(output_unit, '(A)') "All tables converted successfully!"
    end if

contains

    !> Parse command line arguments
    subroutine parse_arguments(input_dir, output_dir, status)
        character(len=*), intent(out) :: input_dir, output_dir
        integer, intent(out) :: status

        integer :: n_args

        n_args = command_argument_count()

        if (n_args == 0) then
            ! Use default directories
            input_dir = "data/tables/cpp_legacy"
            output_dir = "data/tables/fortran_native"
            status = 0
        else if (n_args == 2) then
            ! Use provided directories
            call get_command_argument(1, input_dir)
            call get_command_argument(2, output_dir)
            status = 0
        else
            status = -1
        end if

    end subroutine parse_arguments


    !> Print usage information
    subroutine print_usage()
        write(output_unit, '(A)') ""
        write(output_unit, '(A)') "Usage:"
        write(output_unit, '(A)') "  fpm run convert_tables"
        write(output_unit, '(A)') "  fpm run convert_tables -- <input_dir> <output_dir>"
        write(output_unit, '(A)') ""
        write(output_unit, '(A)') "Default directories:"
        write(output_unit, '(A)') "  Input:  data/tables/cpp_legacy"
        write(output_unit, '(A)') "  Output: data/tables/fortran_native"
        write(output_unit, '(A)') ""
    end subroutine print_usage


    !> Construct C++ table filename
    subroutine construct_cpp_filename(dir, U, filename)
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: U
        character(len=*), intent(out) :: filename

        character(len=32) :: u_string

        ! Format U with proper leading zero
        if (U < 1.0_dp) then
            write(u_string, '(F4.2)') U  ! "0.90"
        else if (U < 10.0_dp) then
            write(u_string, '(F4.2)') U  ! "4.00"
        else
            write(u_string, '(F5.2)') U  ! "10.00"
        end if

        ! Remove leading spaces
        u_string = adjustl(u_string)

        ! Construct full path
        filename = trim(dir) // "/lsda_hub_u" // trim(u_string)

    end subroutine construct_cpp_filename


    !> Construct Fortran binary table filename
    subroutine construct_fortran_filename(dir, U, filename)
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: U
        character(len=*), intent(out) :: filename

        character(len=32) :: u_string

        ! Format U with proper leading zero
        if (U < 1.0_dp) then
            write(u_string, '(F4.2)') U  ! "0.90"
        else if (U < 10.0_dp) then
            write(u_string, '(F4.2)') U  ! "4.00"
        else
            write(u_string, '(F5.2)') U  ! "10.00"
        end if

        ! Remove leading spaces
        u_string = adjustl(u_string)

        ! Construct full path
        filename = trim(dir) // "/xc_table_u" // trim(u_string) // ".dat"

    end subroutine construct_fortran_filename

end program convert_tables