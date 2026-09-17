!> Unit tests for bethe_tables module
program test_bethe_tables
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
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
            test("compute_E_xc_half_filling_u4", test_compute_E_xc_half_filling_u4), &
            test("compute_V_xc_half_filling_u4", test_compute_V_xc_half_filling_u4), &
            test("compute_E_xc_above_half_filling_u4", test_compute_E_xc_above_half_filling_u4), &
            test("compute_E_xc_spin_exchange_u4", test_compute_E_xc_spin_exchange_u4), &
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

    !> Regression test for the finite-size U=4 half-filled Lieb-Wu solution.
    !!
    !! The reference is the C++ U=4 table value after the Hartree contribution
    !! has been removed.  It fails if charge scattering uses k instead of sin(k),
    !! if only one spin population is used for charge roots, or if Hartree remains.
    subroutine test_compute_E_xc_half_filling_u4()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: compute_E_xc, generate_xc_table, grid_params_t
        use lsda_constants, only: dp
        use table_io, only: xc_table_t

        real(dp), parameter :: CPP_EXC = -0.300489823_dp
        real(dp) :: E_xc
        integer :: status
        type(grid_params_t) :: params
        type(xc_table_t) :: table

        E_xc = compute_E_xc(0.5_dp, 0.5_dp, 4.0_dp, 100)
        call check(abs(E_xc - CPP_EXC) < 1.0e-4_dp, &
                   "U=4 half-filled E_xc must match the C++ reference table")

        params%n_min = 1.0_dp
        params%n_max = 1.0_dp
        params%n_points = 1
        params%m_points = 1
        params%L = 100
        call generate_xc_table(4.0_dp, params, table, status)
        call check(status == 0, "U=4 reference table point should generate")
        call check(abs(table%exc(1, 1) - CPP_EXC) < 1.0e-4_dp, &
                   "Generated U=4 table must match the C++ reference table")

        deallocate(table%n_grid, table%m_grid, table%exc, table%vxc_up, table%vxc_down)
    end subroutine test_compute_E_xc_half_filling_u4

    !> Regression test for finite and symmetric U=4 potentials at half filling.
    !!
    !! The discrete derivative must use the one-sided lower-density derivative
    !! at n=1, rather than attempting a Newton solve above half filling.
    subroutine test_compute_V_xc_half_filling_u4()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: compute_V_xc_numerical, xc_potentials_t
        use lsda_constants, only: dp

        type(xc_potentials_t) :: v_xc

        v_xc = compute_V_xc_numerical(0.5_dp, 0.5_dp, 4.0_dp, 100)

        call check(.not. ieee_is_nan(v_xc%v_xc_up), &
                   "U=4 half-filled V_xc_up must be finite")
        call check(.not. ieee_is_nan(v_xc%v_xc_down), &
                   "U=4 half-filled V_xc_down must be finite")
        call check(abs(v_xc%v_xc_up - v_xc%v_xc_down) < 1.0e-10_dp, &
                   "U=4 half-filled V_xc must preserve spin symmetry")
    end subroutine test_compute_V_xc_half_filling_u4

    !> Regression test for particle-hole evaluation above half filling.
    !!
    !! The energies and potentials must remain finite when a finite-difference
    !! derivative evaluates a density above the n=1 table boundary.
    subroutine test_compute_E_xc_above_half_filling_u4()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: compute_E_xc, compute_V_xc_numerical, xc_potentials_t
        use lsda_constants, only: dp

        real(dp) :: e_xc
        type(xc_potentials_t) :: v_xc

        e_xc = compute_E_xc(0.75_dp, 0.75_dp, 4.0_dp, 20)
        v_xc = compute_V_xc_numerical(0.75_dp, 0.75_dp, 4.0_dp, 20)

        call check(.not. ieee_is_nan(e_xc), "Above-half-filled U=4 E_xc must be finite")
        call check(.not. ieee_is_nan(v_xc%v_xc_up), &
                   "Above-half-filled U=4 V_xc_up must be finite")
        call check(.not. ieee_is_nan(v_xc%v_xc_down), &
                   "Above-half-filled U=4 V_xc_down must be finite")
    end subroutine test_compute_E_xc_above_half_filling_u4

    !> Regression test for the minority-spin rapidity convention.
    !!
    !! The negative-m half of a table swaps which physical spin is the
    !! minority.  Both labelings must solve the same Lieb-Wu state and return
    !! a finite, spin-exchange-symmetric energy.
    subroutine test_compute_E_xc_spin_exchange_u4()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: compute_E_xc
        use lsda_constants, only: dp

        real(dp) :: e_xc_up, e_xc_down

        e_xc_up = compute_E_xc(0.1_dp, 0.7_dp, 4.0_dp, 20)
        e_xc_down = compute_E_xc(0.7_dp, 0.1_dp, 4.0_dp, 20)

        call check(.not. ieee_is_nan(e_xc_up), &
                   "Negative-m state must not enter an invalid Bethe sector")
        call check(.not. ieee_is_nan(e_xc_down), &
                   "Positive-m state must produce a finite XC energy")
        call check(abs(e_xc_up - e_xc_down) < 1.0e-12_dp, &
                   "XC energy must be invariant under spin exchange")
    end subroutine test_compute_E_xc_spin_exchange_u4

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
        call check(ABS(params%n_max - 1.0_dp) < TOL, "Default n_max should be 1.0")
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
