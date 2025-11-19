
program test_table_io
    use fortuno_serial, only: execute_serial_cmd_app
    implicit none
    
    call execute_serial_cmd_app(get_table_io_tests())
    
contains

    function get_table_io_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests
        
        tests = test_list([ &
            test("read_cpp_u4", test_read_cpp_u4), &
            test("extract_u_standard", test_extract_u_standard), &
            test("extract_u_binary", test_extract_u_binary), &
            test("extract_u_invalid", test_extract_u_invalid), &
            test("binary_roundtrip", test_binary_roundtrip), &
            test("table_dimensions", test_table_dimensions), &
            test("table_ranges", test_table_ranges), &
            test("missing_file", test_missing_file), &
            test("deallocate_table", test_deallocate_table), &
            test("multiple_tables", test_multiple_tables) &
        ])
    end function get_table_io_tests


    subroutine test_read_cpp_u4()
        use fortuno_serial, only: check => serial_check
        use table_io
        use lsda_constants, only: dp
        
        type(xc_table_t) :: table
        integer :: status

        call read_cpp_table("data/tables/cpp_legacy/lsda_hub_u4.00", table, status)

        call check(status == 0, "Read status should be 0")
        call check(abs(table%U - 4.0_dp) < 1.0e-10_dp, "U should be 4.0")
        call check(table%n_points_n == 50, "Should have 50 density points")
        call check(table%n_points_m > 0, "Should have magnetization points")
        call check(allocated(table%n_grid), "n_grid should be allocated")
        call check(allocated(table%m_grid), "m_grid should be allocated")
        call check(allocated(table%exc), "exc should be allocated")
        call check(allocated(table%vxc_up), "vxc_up should be allocated")
        call check(allocated(table%vxc_down), "vxc_down should be allocated")

        call deallocate_table(table)
    end subroutine test_read_cpp_u4


    subroutine test_extract_u_standard()
        use fortuno_serial, only: check => serial_check
        use table_io, only: extract_U_from_filename
        use lsda_constants, only: dp
        
        real(dp) :: U
        integer :: status

        call extract_U_from_filename("data/tables/cpp_legacy/lsda_hub_u4.00", U, status)
        
        call check(status == 0, "Should extract U from standard filename")
        call check(abs(U - 4.0_dp) < 1.0e-10_dp, "U should be 4.0")
    end subroutine test_extract_u_standard


    subroutine test_extract_u_binary()
        use fortuno_serial, only: check => serial_check
        use table_io, only: extract_U_from_filename
        use lsda_constants, only: dp
        
        real(dp) :: U
        integer :: status

        call extract_U_from_filename("xc_table_u10.00.dat", U, status)
        print *, "DEBUG TEST: status = ", status
        print *, "DEBUG TEST: U = ", U
        
        call check(status == 0, "Should extract U from binary filename")
        call check(abs(U - 10.0_dp) < 1.0e-10_dp, "U should be 10.0")
    end subroutine test_extract_u_binary


    subroutine test_extract_u_invalid()
        use fortuno_serial, only: check => serial_check
        use table_io, only: extract_U_from_filename
        use lsda_constants, only: dp
        
        real(dp) :: U
        integer :: status

        call extract_U_from_filename("invalid_filename.dat", U, status)
        
        call check(status /= 0, "Should fail for invalid filename")
    end subroutine test_extract_u_invalid


    subroutine test_binary_roundtrip()
        use fortuno_serial, only: check => serial_check
        use table_io
        use lsda_constants, only: dp
        
        type(xc_table_t) :: table_in, table_out
        integer :: status
        character(len=256) :: temp_file

        call read_cpp_table("data/tables/cpp_legacy/lsda_hub_u4.00", table_in, status)
        call check(status == 0, "Should read C++ table")

        temp_file = "data/tables/fortran_native/test_temp_u4.00.dat"
        call write_fortran_table(temp_file, table_in, status)
        call check(status == 0, "Should write binary table")

        call read_fortran_table(temp_file, table_out, status)
        call check(status == 0, "Should read binary table")

        call check(abs(table_in%U - table_out%U) < 1.0e-14_dp, &
                   "U should match after round-trip")

        call check(table_in%n_points_n == table_out%n_points_n, &
                   "n_points_n should match")
        call check(table_in%n_points_m == table_out%n_points_m, &
                   "n_points_m should match")

        call check(abs(table_in%n_grid(1) - table_out%n_grid(1)) < 1.0e-14_dp, &
                   "First n value should match")
        call check(abs(table_in%n_grid(50) - table_out%n_grid(50)) < 1.0e-14_dp, &
                   "Last n value should match")
        
        call check(abs(table_in%exc(1,1) - table_out%exc(1,1)) < 1.0e-14_dp, &
                   "First exc value should match")
        call check(abs(table_in%vxc_up(100,25) - table_out%vxc_up(100,25)) < 1.0e-14_dp, &
                   "Middle vxc_up value should match")

        call deallocate_table(table_in)
        call deallocate_table(table_out)
    end subroutine test_binary_roundtrip


    subroutine test_table_dimensions()
        use fortuno_serial, only: check => serial_check
        use table_io
        use lsda_constants, only: dp
        
        type(xc_table_t) :: table
        integer :: status

        call read_cpp_table("data/tables/cpp_legacy/lsda_hub_u4.00", table, status)
        call check(status == 0, "Should read table")

        call check(size(table%n_grid) == table%n_points_n, &
                   "n_grid size should match n_points_n")
        call check(size(table%m_grid, 1) == table%n_points_m, &
                   "m_grid first dimension should match n_points_m")
        call check(size(table%m_grid, 2) == table%n_points_n, &
                   "m_grid second dimension should match n_points_n")
        call check(size(table%exc, 1) == table%n_points_m, &
                   "exc first dimension should match n_points_m")
        call check(size(table%exc, 2) == table%n_points_n, &
                   "exc second dimension should match n_points_n")

        call deallocate_table(table)
    end subroutine test_table_dimensions


    subroutine test_table_ranges()
        use fortuno_serial, only: check => serial_check
        use table_io
        use lsda_constants, only: dp
        
        type(xc_table_t) :: table
        integer :: status
        real(dp) :: n_min, n_max, m_min, m_max

        call read_cpp_table("data/tables/cpp_legacy/lsda_hub_u4.00", table, status)
        call check(status == 0, "Should read table")

        n_min = minval(table%n_grid)
        n_max = maxval(table%n_grid)
        
        call check(n_min > 0.0_dp, "Minimum density should be positive")
        call check(n_max <= 1.0_dp + 1.0e-6_dp, "Maximum density should be <= 1")
        call check(n_min < n_max, "Density should be increasing")

        m_min = minval(table%m_grid)
        m_max = maxval(table%m_grid)
        
        call check(m_min >= 0.0_dp - 1.0e-10_dp, "Minimum magnetization should be >= 0")
        call check(m_max <= n_max + 1.0e-6_dp, "Maximum magnetization should be <= n_max")

        call check(maxval(table%exc) <= 0.0_dp + 1.0e-10_dp, &
                   "XC energy should be negative or zero")

        call check(abs(table%vxc_up(1,1) - table%vxc_down(1,1)) < 1.0e-10_dp, &
                   "At m=0, vxc_up should equal vxc_down")

        call deallocate_table(table)
    end subroutine test_table_ranges


    subroutine test_missing_file()
        use fortuno_serial, only: check => serial_check
        use table_io
        use lsda_constants, only: dp
        
        type(xc_table_t) :: table
        integer :: status

        call read_cpp_table("data/tables/cpp_legacy/nonexistent_file.dat", table, status)
        
        call check(status /= 0, "Should fail for missing file")
        call check(.not. allocated(table%n_grid), "Arrays should not be allocated on error")
    end subroutine test_missing_file


    subroutine test_deallocate_table()
        use fortuno_serial, only: check => serial_check
        use table_io
        use lsda_constants, only: dp
        
        type(xc_table_t) :: table
        integer :: status

        call read_cpp_table("data/tables/cpp_legacy/lsda_hub_u4.00", table, status)
        call check(status == 0, "Should read table")
        call check(allocated(table%n_grid), "Arrays should be allocated")

        call deallocate_table(table)
        
        call check(.not. allocated(table%n_grid), "n_grid should be deallocated")
        call check(.not. allocated(table%m_grid), "m_grid should be deallocated")
        call check(.not. allocated(table%exc), "exc should be deallocated")
        call check(.not. allocated(table%vxc_up), "vxc_up should be deallocated")
        call check(.not. allocated(table%vxc_down), "vxc_down should be deallocated")

        call deallocate_table(table)
    end subroutine test_deallocate_table


    subroutine test_multiple_tables()
        use fortuno_serial, only: check => serial_check
        use table_io
        use lsda_constants, only: dp
        
        type(xc_table_t) :: table1, table2, table3
        integer :: status

        call read_cpp_table("data/tables/cpp_legacy/lsda_hub_u2.00", table1, status)
        call check(status == 0, "Should read U=2 table")
        
        call read_cpp_table("data/tables/cpp_legacy/lsda_hub_u4.00", table2, status)
        call check(status == 0, "Should read U=4 table")
        
        call read_cpp_table("data/tables/cpp_legacy/lsda_hub_u10.00", table3, status)
        call check(status == 0, "Should read U=10 table")

        call check(abs(table1%U - 2.0_dp) < 1.0e-10_dp, "table1 U should be 2.0")
        call check(abs(table2%U - 4.0_dp) < 1.0e-10_dp, "table2 U should be 4.0")
        call check(abs(table3%U - 10.0_dp) < 1.0e-10_dp, "table3 U should be 10.0")

        call check(table1%n_points_n > 0, "table1 should have density points")
        call check(table2%n_points_n > 0, "table2 should have density points")
        call check(table3%n_points_n > 0, "table3 should have density points")
        
        call check(table1%n_points_m > 0, "table1 should have magnetization points")
        call check(table2%n_points_m > 0, "table2 should have magnetization points")
        call check(table3%n_points_m > 0, "table3 should have magnetization points")

        call deallocate_table(table1)
        call deallocate_table(table2)
        call deallocate_table(table3)
    end subroutine test_multiple_tables

end program test_table_io