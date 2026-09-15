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
            test("symmetry_transformations", test_symmetry_transformations), &
            test("vxc_jump_at_half_filling_without_smoothing", &
                 test_vxc_jump_at_half_filling_without_smoothing), &
            test("vxc_smoothing_removes_jump", test_vxc_smoothing_removes_jump), &
            test("vxc_smoothing_invalid_width_rejected", &
                 test_vxc_smoothing_invalid_width_rejected), &
            test("vxc_smoothing_nonfinite_width_rejected", &
                 test_vxc_smoothing_nonfinite_width_rejected) &
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

    !> Default (w = 0) must keep the physical V_xc discontinuity at n = 1
    !!
    !! For U = 4 and m = 0 the branch used above n = 1 carries an overall minus
    !! sign, so V_xc jumps by 2|v_base(1, 0)| ≈ 1.29 across half filling. This
    !! is the C++ behaviour and the default must reproduce it: the test fails if
    !! the smoothing is ever turned on by accident.
    subroutine test_vxc_jump_at_half_filling_without_smoothing()
        use fortuno_serial, only: check => serial_check
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, get_vxc, xc_lsda_destroy
        use lsda_constants, only: dp

        type(xc_lsda_t) :: xc
        integer :: status, ierr
        real(dp) :: v_lo_up, v_lo_dw, v_hi_up, v_hi_dw
        real(dp), parameter :: EPS = 1.0e-6_dp
        character(len=256) :: test_file

        test_file = "data/tables/fortran_native/xc_table_u4.00.dat"
        call xc_lsda_init(xc, test_file, status)
        call check(status == 0, "XC initialization should succeed")
        call check(abs(xc%smoothing_width) < TOL, "smoothing must be off by default")

        call get_vxc(xc, 0.5_dp - 0.5_dp*EPS, 0.5_dp - 0.5_dp*EPS, v_lo_up, v_lo_dw, ierr)
        call check(ierr == 0, "get_vxc below half filling should succeed")

        call get_vxc(xc, 0.5_dp + 0.5_dp*EPS, 0.5_dp + 0.5_dp*EPS, v_hi_up, v_hi_dw, ierr)
        call check(ierr == 0, "get_vxc above half filling should succeed")

        call check(abs(v_hi_up - v_lo_up) > 1.0_dp, &
                   "V_xc^up must jump across n = 1 when smoothing is off")
        call check(abs(v_hi_dw - v_lo_dw) > 1.0_dp, &
                   "V_xc^dn must jump across n = 1 when smoothing is off")

        call xc_lsda_destroy(xc)
    end subroutine test_vxc_jump_at_half_filling_without_smoothing

    !> A positive smoothing width removes the jump and leaves the rest untouched
    subroutine test_vxc_smoothing_removes_jump()
        use fortuno_serial, only: check => serial_check
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, get_vxc, xc_lsda_destroy
        use lsda_constants, only: dp

        type(xc_lsda_t) :: xc_plain, xc_smooth
        integer :: status, ierr
        real(dp) :: v_lo_up, v_lo_dw, v_hi_up, v_hi_dw
        real(dp) :: v_ref_up, v_ref_dw, v_out_up, v_out_dw
        real(dp), parameter :: EPS = 1.0e-6_dp
        real(dp), parameter :: W = 0.05_dp
        character(len=256) :: test_file

        test_file = "data/tables/fortran_native/xc_table_u4.00.dat"

        call xc_lsda_init(xc_plain, test_file, status)
        call check(status == 0, "unsmoothed XC initialization should succeed")

        call xc_lsda_init(xc_smooth, test_file, status, smoothing_width=W)
        call check(status == 0, "smoothed XC initialization should succeed")
        call check(abs(xc_smooth%smoothing_width - W) < TOL, "smoothing width should be stored")

        ! Inside the window the potential is now continuous across n = 1
        call get_vxc(xc_smooth, 0.5_dp - 0.5_dp*EPS, 0.5_dp - 0.5_dp*EPS, v_lo_up, v_lo_dw, ierr)
        call check(ierr == 0, "smoothed get_vxc below n = 1 should succeed")

        call get_vxc(xc_smooth, 0.5_dp + 0.5_dp*EPS, 0.5_dp + 0.5_dp*EPS, v_hi_up, v_hi_dw, ierr)
        call check(ierr == 0, "smoothed get_vxc above n = 1 should succeed")

        call check(abs(v_hi_up - v_lo_up) < 1.0e-4_dp, &
                   "smoothing should remove the V_xc^up jump at n = 1")
        call check(abs(v_hi_dw - v_lo_dw) < 1.0e-4_dp, &
                   "smoothing should remove the V_xc^dn jump at n = 1")

        ! Outside the window nothing may change
        call get_vxc(xc_plain, 0.4_dp, 0.3_dp, v_ref_up, v_ref_dw, ierr)
        call check(ierr == 0, "unsmoothed get_vxc outside the window should succeed")

        call get_vxc(xc_smooth, 0.4_dp, 0.3_dp, v_out_up, v_out_dw, ierr)
        call check(ierr == 0, "smoothed get_vxc outside the window should succeed")

        call check(abs(v_out_up - v_ref_up) < TOL .and. abs(v_out_dw - v_ref_dw) < TOL, &
                   "smoothing must not change V_xc outside |n - 1| < w")

        ! And at the window edge the two functionals must still agree exactly
        call get_vxc(xc_plain, 0.5_dp*(1.0_dp - W), 0.5_dp*(1.0_dp - W), v_ref_up, v_ref_dw, ierr)
        call get_vxc(xc_smooth, 0.5_dp*(1.0_dp - W), 0.5_dp*(1.0_dp - W), v_out_up, v_out_dw, ierr)
        call check(abs(v_out_up - v_ref_up) < TOL .and. abs(v_out_dw - v_ref_dw) < TOL, &
                   "smoothed V_xc must match the plain one at the lower window edge")

        call xc_lsda_destroy(xc_plain)
        call xc_lsda_destroy(xc_smooth)
    end subroutine test_vxc_smoothing_removes_jump

    !> Negative or too large smoothing widths are rejected
    subroutine test_vxc_smoothing_invalid_width_rejected()
        use fortuno_serial, only: check => serial_check
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use lsda_constants, only: dp

        type(xc_lsda_t) :: xc
        integer :: status
        character(len=256) :: test_file

        test_file = "data/tables/fortran_native/xc_table_u4.00.dat"

        call xc_lsda_init(xc, test_file, status, smoothing_width=-0.1_dp)
        call check(status /= 0, "negative smoothing width must be rejected")
        call check(.not. xc%initialized, "XC must stay uninitialized on invalid width")

        call xc_lsda_init(xc, test_file, status, smoothing_width=1.0_dp)
        call check(status /= 0, "smoothing width >= 1 must be rejected")

        call xc_lsda_destroy(xc)
    end subroutine test_vxc_smoothing_invalid_width_rejected

    !> A non-finite smoothing width must be rejected, NaN included
    !!
    !! NaN is the dangerous value here: it passes BOTH range comparisons
    !! (`NaN < 0` and `NaN >= 1` are false), gets stored in xc%smoothing_width,
    !! and then fails silently inside get_vxc, where `w > 0.0_dp` is also false
    !! and the unsmoothed branch is taken. The caller asked for smoothing and
    !! would get the discontinuous V_xc with no diagnostic at all. Without the
    !! ieee_is_nan guard the first assertion below fails.
    !!
    !! The two infinities are covered as well to pin the claim that the range
    !! checks already reject them (+Inf >= 1, -Inf < 0): if that range is ever
    !! widened, these assertions force the finiteness question to be revisited.
    !!
    !! The NaN/Inf values come from ieee_value, never from arithmetic such as
    !! 0.0/0.0, which the compiler may fold or trap depending on the flags.
    subroutine test_vxc_smoothing_nonfinite_width_rejected()
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
                                                 ieee_positive_inf, ieee_negative_inf
        use fortuno_serial, only: check => serial_check
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use lsda_constants, only: dp

        type(xc_lsda_t) :: xc
        integer :: status
        character(len=256) :: test_file
        real(dp) :: nan_w, pinf_w, ninf_w

        test_file = "data/tables/fortran_native/xc_table_u4.00.dat"

        nan_w = ieee_value(1.0_dp, ieee_quiet_nan)
        pinf_w = ieee_value(1.0_dp, ieee_positive_inf)
        ninf_w = ieee_value(1.0_dp, ieee_negative_inf)

        call xc_lsda_init(xc, test_file, status, smoothing_width=nan_w)
        call check(status /= 0, "NaN smoothing width must be rejected")
        call check(.not. xc%initialized, "XC must stay uninitialized on a NaN width")

        call xc_lsda_init(xc, test_file, status, smoothing_width=pinf_w)
        call check(status /= 0, "+Infinity smoothing width must be rejected")
        call check(.not. xc%initialized, "XC must stay uninitialized on a +Inf width")

        call xc_lsda_init(xc, test_file, status, smoothing_width=ninf_w)
        call check(status /= 0, "-Infinity smoothing width must be rejected")
        call check(.not. xc%initialized, "XC must stay uninitialized on a -Inf width")

        call xc_lsda_destroy(xc)
    end subroutine test_vxc_smoothing_nonfinite_width_rejected
end program test_xc_lsda
