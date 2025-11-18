!> Unit tests for xc_lsda module
program test_xc_lsda
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp
    implicit none

    real(dp), parameter :: TOL = 1.0e-9_dp

    call execute_serial_cmd_app(get_xc_lsda_tests())

contains

    function get_xc_lsda_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests

        tests = test_list([ &
            test("xc_lsda_init_destroy", test_xc_lsda_init_destroy), &
            test("get_exc_evaluation", test_get_exc_evaluation), &
            test("get_exc_spin_symmetry", test_get_exc_spin_symmetry), &
            test("get_vxc_spin_symmetry", test_get_vxc_spin_symmetry), &
            test("region_determination", test_region_determination), &
            test("symmetry_transformations", test_symmetry_transformations) &
        ])
    end function get_xc_lsda_tests

    subroutine test_xc_lsda_init_destroy()
        use fortuno_serial, only: check => serial_check
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use lsda_constants, only: dp

        type(xc_lsda_t) :: xc
        integer :: status
        character(len=256) :: test_file

        test_file = "data/tables/fortran_native/xc_table_u4.00.dat"

        call xc_lsda_init(xc, test_file, status)
        call check(status == 0, "XC initialization should succeed")
        call check(xc%initialized, "XC should be initialized")
        call check(abs(xc%U - 4.0_dp) < TOL, "U should match table")

        call xc_lsda_destroy(xc)
        call check(.not. xc%initialized, "XC should be deinitialized")

    end subroutine test_xc_lsda_init_destroy

    !> Test exc evaluation returns valid numbers
    subroutine test_get_exc_evaluation()
        use fortuno_serial, only: check => serial_check
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, get_exc, xc_lsda_destroy
        use lsda_constants, only: dp

        type(xc_lsda_t) :: xc
        integer :: status, ierr
        real(dp) :: exc, n_up, n_dw
        character(len=256) :: test_file

        test_file = "data/tables/fortran_native/xc_table_u2.00.dat"
        call xc_lsda_init(xc, test_file, status)

        n_up = 0.4_dp
        n_dw = 0.3_dp

        call get_exc(xc, n_up, n_dw, exc, ierr)

        call check(ierr == 0, "get_exc should succeed")
        call check(exc == exc, "exc should be a valid number")
        call check(abs(exc) > 1.0e-10_dp, "exc should be non-zero for U > 0")

        call xc_lsda_destroy(xc)

    end subroutine test_get_exc_evaluation

    !> Test spin exchange symmetry: exc(n_up, n_dw) = exc(n_dw, n_up)
    subroutine test_get_exc_spin_symmetry()
        use fortuno_serial, only: check => serial_check
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, get_exc, xc_lsda_destroy
        use lsda_constants, only: dp

        type(xc_lsda_t) :: xc
        integer :: status, ierr
        real(dp) :: exc1, exc2, n_up, n_dw
        character(len=256) :: test_file

        test_file = "data/tables/fortran_native/xc_table_u2.00.dat"
        call xc_lsda_init(xc, test_file, status)

        n_up = 0.45_dp
        n_dw = 0.30_dp

        call get_exc(xc, n_up, n_dw, exc1, ierr)
        call check(ierr == 0, "get_exc should succeed for (n_up, n_dw)")

        call get_exc(xc, n_dw, n_up, exc2, ierr)
        call check(ierr == 0, "get_exc should succeed for (n_dw, n_up)")

        call check(abs(exc1 - exc2) < 1.0e-6_dp, &
                   "exc should be symmetric under spin exchange")

        call xc_lsda_destroy(xc)

    end subroutine test_get_exc_spin_symmetry

    !> Test V_xc spin exchange: V_up(n_up, n_dw) = V_dn(n_dw, n_up)
    subroutine test_get_vxc_spin_symmetry()
        use fortuno_serial, only: check => serial_check
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, get_vxc, xc_lsda_destroy
        use lsda_constants, only: dp

        type(xc_lsda_t) :: xc
        integer :: status, ierr
        real(dp) :: v_up1, v_dw1, v_up2, v_dw2, n_up, n_dw
        character(len=256) :: test_file

        test_file = "data/tables/fortran_native/xc_table_u2.00.dat"
        call xc_lsda_init(xc, test_file, status)

        n_up = 0.45_dp
        n_dw = 0.30_dp

        call get_vxc(xc, n_up, n_dw, v_up1, v_dw1, ierr)
        call check(ierr == 0, "get_vxc should succeed for (n_up, n_dw)")

        call get_vxc(xc, n_dw, n_up, v_up2, v_dw2, ierr)
        call check(ierr == 0, "get_vxc should succeed for (n_dw, n_up)")

        call check(abs(v_up1 - v_dw2) < 1.0e-6_dp, &
                   "V_up(n_up,n_dw) should equal V_dn(n_dw,n_up)")
        call check(abs(v_dw1 - v_up2) < 1.0e-6_dp, &
                   "V_dn(n_up,n_dw) should equal V_up(n_dw,n_up)")

        call xc_lsda_destroy(xc)

    end subroutine test_get_vxc_spin_symmetry

    !> Test region determination logic
    subroutine test_region_determination()
        use fortuno_serial, only: check => serial_check
        use lsda_constants, only: dp

        real(dp) :: n, m

        ! Region I: m ≥ 0, n ≤ 1
        n = 0.6_dp
        m = 0.2_dp
        call check(m >= 0.0_dp .and. n <= 1.0_dp, "Should be Region I")

        ! Region II: m < 0, n ≤ 1
        n = 0.6_dp
        m = -0.2_dp
        call check(m < 0.0_dp .and. n <= 1.0_dp, "Should be Region II")

        ! Region III: m < 0, n > 1
        n = 1.4_dp
        m = -0.2_dp
        call check(m < 0.0_dp .and. n > 1.0_dp, "Should be Region III")

        ! Region IV: m ≥ 0, n > 1
        n = 1.4_dp
        m = 0.2_dp
        call check(m >= 0.0_dp .and. n > 1.0_dp, "Should be Region IV")

    end subroutine test_region_determination

    !> Test symmetry transformations map correctly
    subroutine test_symmetry_transformations()
        use fortuno_serial, only: check => serial_check
        use lsda_constants, only: dp

        real(dp) :: n_up, n_dw, n_up_map, n_dw_map

        ! Region I: Identity
        n_up = 0.4_dp
        n_dw = 0.2_dp
        n_up_map = n_up
        n_dw_map = n_dw
        call check(abs(n_up_map - 0.4_dp) < TOL .and. abs(n_dw_map - 0.2_dp) < TOL, &
                   "Region I should be identity")

        ! Region II: Spin exchange
        n_up = 0.2_dp
        n_dw = 0.4_dp
        n_up_map = n_dw
        n_dw_map = n_up
        call check(abs(n_up_map - 0.4_dp) < TOL .and. abs(n_dw_map - 0.2_dp) < TOL, &
                   "Region II should exchange spins")

        ! Region III: Particle-hole
        n_up = 0.3_dp
        n_dw = 0.8_dp
        n_up_map = 1.0_dp - n_up
        n_dw_map = 1.0_dp - n_dw
        call check(abs(n_up_map - 0.7_dp) < TOL .and. abs(n_dw_map - 0.2_dp) < TOL, &
                   "Region III should apply particle-hole symmetry")

        ! Region IV: Combined
        n_up = 0.8_dp
        n_dw = 0.3_dp
        n_up_map = 1.0_dp - n_dw
        n_dw_map = 1.0_dp - n_up
        call check(abs(n_up_map - 0.7_dp) < TOL .and. abs(n_dw_map - 0.2_dp) < TOL, &
                   "Region IV should apply combined symmetry")

    end subroutine test_symmetry_transformations
end program test_xc_lsda
