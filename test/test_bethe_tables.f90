!> Unit tests for bethe_tables module
program test_bethe_tables
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp
    implicit none

    real(dp), parameter :: TOL = 1.0e-9_dp
    
    call execute_serial_cmd_app(get_bethe_tables_tests())
    
contains

    function get_bethe_tables_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests
        
        tests = test_list([ &
            test("compute_E0_half_filling", test_compute_E0_half_filling), &
            test("compute_E0_polarized", test_compute_E0_polarized), &
            test("compute_E_xc_U0", test_compute_E_xc_U0), &
            test("compute_V_xc_symmetric", test_compute_V_xc_symmetric), &
            test("grid_params_defaults", test_grid_params_defaults), &
            test("generate_small_table", test_generate_small_table) &
        ])
    end function get_bethe_tables_tests

    !> Test E0 at half-filling (n_up = n_dn = 0.5)
    subroutine test_compute_E0_half_filling()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: compute_E0
        use lsda_constants, only: dp
        
        real(dp) :: E0, n_up, n_dw
        integer :: L
        
        n_up = 0.5_dp
        n_dw = 0.5_dp
        L = 100
        
        E0 = compute_E0(n_up, n_dw, L)
        
        ! At half-filling with N=100, E0 should be negative
        call check(E0 < 0.0_dp, "E0 should be negative at half-filling")
        
        ! Check reasonable magnitude (kinetic energy ~ -4t for N particles)
        call check(abs(E0) < 400.0_dp, "E0 magnitude should be reasonable")
        
    end subroutine test_compute_E0_half_filling

    !> Test E0 for fully polarized case (n_dw = 0)
    subroutine test_compute_E0_polarized()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: compute_E0
        use lsda_constants, only: dp
        
        real(dp) :: E0, n_up, n_dw
        integer :: L
        
        n_up = 0.5_dp
        n_dw = 0.0_dp
        L = 100
        
        E0 = compute_E0(n_up, n_dw, L)
        
        ! Should still be negative
        call check(E0 < 0.0_dp, "E0 should be negative for polarized case")
        
    end subroutine test_compute_E0_polarized

    !> Test that E_xc = 0 for U = 0 (within tolerance)
    subroutine test_compute_E_xc_U0()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: compute_E_xc
        use lsda_constants, only: dp
        
        real(dp) :: E_xc, n_up, n_dw, U
        integer :: L
        
        n_up = 0.3_dp
        n_dw = 0.2_dp
        U = 0.0_dp
        L = 50
        
        E_xc = compute_E_xc(n_up, n_dw, U, L)
        print *, "E_xc for U=0:", E_xc
        
        ! For U=0, E_xc should be zero (E_BA = E0)
        call check(abs(E_xc) < 1.0e-6_dp, "E_xc should be ~0 for U=0")
        
    end subroutine test_compute_E_xc_U0

    !> Test V_xc symmetry: for n_up = n_dw, V_xc_up = V_xc_dw
    subroutine test_compute_V_xc_symmetric()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: compute_V_xc_numerical, xc_potentials_t
        use lsda_constants, only: dp
        
        real(dp) :: n_up, n_dw, U
        integer :: L
        type(xc_potentials_t) :: v_xc
        
        n_up = 0.25_dp
        n_dw = 0.25_dp  ! Symmetric
        U = 1.0_dp
        L = 20

        v_xc = compute_V_xc_numerical(n_up, n_dw, U, L)

        ! For symmetric case, potentials should be equal
        call check(abs(v_xc%v_xc_up - v_xc%v_xc_down) < 1.0e-6_dp, &
                   "V_xc_up should equal V_xc_down for symmetric case")
        
    end subroutine test_compute_V_xc_symmetric

    !> Test grid_params_t default values
    subroutine test_grid_params_defaults()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: grid_params_t
        use lsda_constants, only: dp
        
        type(grid_params_t) :: params
        
        ! Check defaults
        call check(ABS(params%n_min - 0.1_dp) < TOL, "Default n_min should be 0.1")
        call check(ABS(params%n_max - 2.0_dp) < TOL, "Default n_max should be 2.0")
        call check(params%n_points == 50, "Default n_points should be 50")
        call check(params%m_points == 51, "Default m_points should be 51")
        call check(params%L == 100, "Default L should be 100")
        
    end subroutine test_grid_params_defaults

    !> Test generating a small XC table
    subroutine test_generate_small_table()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: generate_xc_table, grid_params_t
        use table_io, only: xc_table_t
        use lsda_constants, only: dp
        
        type(grid_params_t) :: params
        type(xc_table_t) :: table
        real(dp) :: U
        integer :: status
        
        ! Small grid for fast testing
        params%n_min = 0.25_dp
        params%n_max = 0.25_dp
        params%n_points = 1
        params%m_points = 1
        params%L = 8
        U = 0.1_dp
        
        call generate_xc_table(U, params, table, status)
        
        call check(status == 0, "Table generation should succeed")
        call check(abs(table%U - U) < 1.0e-10_dp, "Table U should match input")
        call check(table%n_points_n == 1, "Should have 3 density points")
        call check(table%n_points_m == 1, "Should have 5 magnetization points")
        call check(allocated(table%exc), "exc array should be allocated")
        call check(allocated(table%vxc_up), "vxc_up array should be allocated")
        call check(allocated(table%vxc_down), "vxc_down array should be allocated")

        ! Clean up
        if (allocated(table%n_grid)) deallocate(table%n_grid)
        if (allocated(table%m_grid)) deallocate(table%m_grid)
        if (allocated(table%exc)) deallocate(table%exc)
        if (allocated(table%vxc_up)) deallocate(table%vxc_up)
        if (allocated(table%vxc_down)) deallocate(table%vxc_down)
        
    end subroutine test_generate_small_table

end program test_bethe_tables