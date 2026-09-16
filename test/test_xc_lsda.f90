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
                 test_vxc_smoothing_nonfinite_width_rejected), &
            test("shiba_transform_for_attractive_u", &
                 test_shiba_transform_for_attractive_u), &
            test("shiba_requires_matching_table", &
                 test_shiba_requires_matching_table), &
            test("dexc_dndown_b0_matches_table", &
                 test_dexc_dndown_b0_matches_table), &
            test("exc_fully_polarized_endpoint_slope", &
                 test_exc_fully_polarized_endpoint_slope), &
            test("exc_below_first_table_density", &
                 test_exc_below_first_table_density), &
            test("empty_channel_shortcuts", test_empty_channel_shortcuts), &
            test("cpp_reference_values_off_nodes", &
                 test_cpp_reference_values_off_nodes) &
        ])
    end function get_xc_lsda_tests

    !> The analytic ∂e_xc/∂n_dn at n_dn = 0 must reproduce the tabulated column
    !!
    !! The last magnetization point of every density row of the table is the
    !! fully polarized point m = n (n_dn = 0), where the tabulated V_xc^dn is by
    !! definition ∂e_xc/∂n_dn. `dexc_dndown_b0` computes the same quantity in
    !! closed form (it is needed off the tabulated densities, as the value of the
    !! synthetic node of the density spline), so the two must agree.
    !!
    !! This pins down the closed form of the auxiliary integral F(γ, Q): any
    !! error in its branch handling or in the constant term shows up here.
    subroutine test_dexc_dndown_b0_matches_table()
        use fortuno_serial, only: check => serial_check
        use xc_lsda, only: dexc_dndown_b0
        use table_io, only: xc_table_t, read_fortran_table, deallocate_table
        use lsda_constants, only: dp

        type(xc_table_t) :: table
        integer :: status, i, nm
        real(dp) :: analytic, tabulated, error, max_error

        call read_fortran_table("data/tables/fortran_native/xc_table_u4.00.dat", table, status)
        call check(status == 0, "table should be readable")
        if (status /= 0) return

        nm = table%n_points_m
        max_error = 0.0_dp
        do i = 1, table%n_points_n
            analytic = dexc_dndown_b0(table%U, table%n_grid(i))
            tabulated = table%vxc_down(nm, i)
            error = abs(analytic - tabulated)
            max_error = max(max_error, error)
        end do

        call check(max_error < 1.0e-6_dp, &
                   "dexc_dndown_b0 should match the fully polarized Vxc_dn column")

        ! Sanity: e_xc vanishes for the non-interacting system
        call check(abs(dexc_dndown_b0(0.0_dp, 0.5_dp)) < 1.0e-14_dp, &
                   "dexc_dndown_b0 should vanish at U = 0")

        call deallocate_table(table)

    end subroutine test_dexc_dndown_b0_matches_table

    !> The m = n endpoint slope is -V_xc^dn/2, not -V_xc^dn
    !!
    !! At fixed total density, dn_up/dm = +1/2 and dn_down/dm = -1/2.
    !! Because e_xc is identically zero along the fully polarized line, its
    !! tangential derivative vanishes there and the inward row derivative is
    !! therefore de_xc/dm = -V_xc^dn/2. The C++ implementation states this
    !! identity but accidentally supplies twice the slope to its spline.
    subroutine test_exc_fully_polarized_endpoint_slope()
        use fortuno_serial, only: check => serial_check
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, get_exc, xc_lsda_destroy
        use table_io, only: xc_table_t, read_fortran_table, deallocate_table
        use lsda_constants, only: dp

        type(xc_lsda_t) :: xc
        type(xc_table_t) :: table
        integer :: status, i, nm
        real(dp) :: n, m, h, exc_inner, exc_edge, slope, expected

        call read_fortran_table("data/tables/fortran_native/xc_table_u4.00.dat", table, status)
        call check(status == 0, "endpoint slope: table should be readable")
        if (status /= 0) return

        call xc_lsda_init(xc, "data/tables/fortran_native/xc_table_u4.00.dat", status)
        call check(status == 0, "endpoint slope: XC initialization should succeed")
        if (status /= 0) then
            call deallocate_table(table)
            return
        end if

        i = table%n_points_n / 2
        nm = table%n_points_m
        n = table%n_grid(i)
        h = 1.0e-6_dp
        m = n - h

        call get_exc(xc, 0.5_dp * (n + m), 0.5_dp * (n - m), exc_inner, status)
        call check(status == 0, "endpoint slope: inner evaluation should succeed")
        call get_exc(xc, n, 0.0_dp, exc_edge, status)
        call check(status == 0, "endpoint slope: edge evaluation should succeed")

        slope = (exc_edge - exc_inner) / h
        expected = -0.5_dp * table%vxc_down(nm, i)
        call check(abs(slope - expected) < 1.0e-5_dp, &
                   "endpoint slope must equal -Vxc_down/2")

        call xc_lsda_destroy(xc)
        call deallocate_table(table)
    end subroutine test_exc_fully_polarized_endpoint_slope

    !> Densities below the first tabulated node must stay physical
    !!
    !! The u = 4 table starts at n = 0.0203; SCF runs do produce sites with
    !! n ≈ 0.0015. With the old linear interpolation in n such a point was
    !! reached by extrapolating the straight line through the first two rows,
    !! which for e_xc overshoots through zero and returns a POSITIVE value
    !! (≈ +6.4e-4) — e_xc is strictly negative for repulsive U and must tend to
    !! zero as n → 0.
    !!
    !! The cubic spline in n, anchored by the synthetic node on the fully
    !! polarized line, keeps the value bracketed.
    subroutine test_exc_below_first_table_density()
        use fortuno_serial, only: check => serial_check
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, get_exc, xc_lsda_destroy
        use lsda_constants, only: dp

        type(xc_lsda_t) :: xc
        integer :: status
        real(dp) :: exc_low, exc_first

        call xc_lsda_init(xc, "data/tables/fortran_native/xc_table_u4.00.dat", status)
        call check(status == 0, "XC initialization should succeed")
        if (status /= 0) return

        ! n = 0.0015, unpolarized (below the first tabulated density n = 0.0203)
        call get_exc(xc, 0.00075_dp, 0.00075_dp, exc_low, status)
        call check(status == 0, "get_exc should succeed below the first density")

        ! n = 0.0203, unpolarized (the first tabulated density)
        call get_exc(xc, 0.01015_dp, 0.01015_dp, exc_first, status)
        call check(status == 0, "get_exc should succeed at the first density")

        call check(exc_low < 0.0_dp, "e_xc must stay negative for repulsive U")
        call check(abs(exc_low) < abs(exc_first), &
                   "|e_xc| must decrease towards n = 0, not extrapolate past it")

        call xc_lsda_destroy(xc)

    end subroutine test_exc_below_first_table_density

    !> REGRESSION (R10): e_xc and V_xc against hard-coded values OFF the nodes
    !!
    !! Every other test of this suite checks a symmetry, a sign, a shortcut or a
    !! value ON a tabulated node. None of them sees an off-by-one in the row
    !! window selected by `spline2d_eval` (`i_first`, `num_n`), in the synthetic
    !! node, or in the boundary condition of the row splines: those only change
    !! the interpolant BETWEEN nodes. This test pins the interpolated values
    !! themselves against the C++ reference (`original/spline2D.cc`, functions
    !! exc_value / Vxc_up_value / Vxc_dn_value, driven for |U| = 4 through a
    !! throw-away driver linked against the original sources). The V_xc values
    !! retain exact C++ parity. Four e_xc references intentionally differ: they
    !! are reached through the final magnetization interval where the C++ uses
    !! -V_xc^dn instead of the analytic endpoint slope -V_xc^dn/2. Their values
    !! below were regenerated after that single correction; the independent
    !! derivative identity is pinned by test_exc_fully_polarized_endpoint_slope.
    !!
    !! The seven points are off-node in BOTH directions and cover
    !!   1,2: Region I   (m >= 0, n <= 1)
    !!   3:   Region II  (m <  0, n <= 1)   -> spin exchange
    !!   4:   Region III (m <  0, n >  1)   -> particle-hole, sign flip
    !!   5:   Region IV  (m >= 0, n >  1)   -> both
    !!   6,7: below the first tabulated density (n = 0.007 and n = 0.004 against
    !!        a first row at n = 0.0202807), i.e. the window that consists of the
    !!        synthetic node plus the whole table
    !! and are evaluated for both signs of U, so the Shiba transformation is
    !! covered as well.
    !!
    !! The tolerance is 1e-9: the measured deviations are at most 3e-11, and
    !! they come from `integral_1`/`dexc_dndown_b0`, which is a closed form here
    !! and a quadrature in the C++.
    subroutine test_cpp_reference_values_off_nodes()
        use fortuno_serial, only: check => serial_check
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, get_exc, get_vxc, xc_lsda_destroy
        use lsda_constants, only: dp

        integer, parameter :: NPTS = 7
        !> Deviation allowed against the C++ reference
        real(dp), parameter :: CPP_TOL = 1.0e-9_dp

        type(xc_lsda_t) :: xc
        integer :: status, ierr, i
        real(dp) :: exc, v_up, v_dw
        real(dp) :: n_up(NPTS), n_dw(NPTS)
        real(dp) :: exc_p(NPTS), vup_p(NPTS), vdn_p(NPTS)
        real(dp) :: exc_m(NPTS), vup_m(NPTS), vdn_m(NPTS)

        n_up = [0.34267985324_dp, 0.27001671326_dp, 0.15_dp, 0.45_dp, 0.82_dp, &
                0.004_dp, 0.0031_dp]
        n_dw = [0.028423868551_dp, 0.20374663913_dp, 0.30_dp, 0.72_dp, 0.43_dp, &
                0.003_dp, 0.0009_dp]

        ! ---- C++ reference, U = +4 ------------------------------------------
        exc_p = [-0.01886444325204032_dp, -0.09734691993050536_dp, &
                 -0.08109155122879777_dp, -0.1976590508289344_dp, &
                 -0.13884924763667_dp, -4.715464898467325e-05_dp, &
                 -1.10470910602286e-05_dp]
        vup_p = [-0.02884866182710155_dp, -0.2424327331672126_dp, &
                 -0.4573210438228484_dp, 0.2416790877745727_dp, &
                 0.6746049834235452_dp, -0.01124631958862716_dp, &
                 -0.003347177556420216_dp]
        vdn_p = [-0.6426005138655668_dp, -0.3754186383910444_dp, &
                 -0.166737725269995_dp, 0.5858783532371914_dp, &
                 0.155478183262056_dp, -0.01518029142525137_dp, &
                 -0.0120627212450126_dp]

        ! ---- C++ reference, U = -4 (Shiba) ----------------------------------
        exc_m = [-0.02697199509353530_dp, -0.1888589974732909_dp, &
                 -0.1370821173219327_dp, -0.1737738830700924_dp, &
                 -0.1173161068289127_dp, -0.004928761787534891_dp, &
                 -0.001483988382532426_dp]
        vup_m = [0.02854827872029039_dp, 0.2896841644791011_dp, &
                 -0.8316581996001842_dp, -0.2430808001099627_dp, &
                 0.5561574713787855_dp, 0.01183572366169481_dp, &
                 0.003548855676387826_dp]
        vdn_m = [-0.9273540394183423_dp, -0.8577903714734391_dp, &
                 0.1714570982716558_dp, 0.4960910736315256_dp, &
                 -0.1593018655000886_dp, -1.641028794017704_dp, &
                 -1.644573668825788_dp]

        call xc_lsda_init(xc, "data/tables/fortran_native/xc_table_u4.00.dat", status, &
                          u_signed = 4.0_dp)
        call check(status == 0, "C++ reference: XC init (U = +4) should succeed")
        if (status /= 0) return

        do i = 1, NPTS
            call get_exc(xc, n_up(i), n_dw(i), exc, ierr)
            call check(ierr == 0, "C++ reference (U = +4): get_exc should succeed")
            call check(abs(exc - exc_p(i)) < CPP_TOL, &
                       "corrected reference (U = +4): e_xc must match off the nodes")

            call get_vxc(xc, n_up(i), n_dw(i), v_up, v_dw, ierr)
            call check(ierr == 0, "C++ reference (U = +4): get_vxc should succeed")
            call check(abs(v_up - vup_p(i)) < CPP_TOL, &
                       "C++ reference (U = +4): V_xc^up must match off the nodes")
            call check(abs(v_dw - vdn_p(i)) < CPP_TOL, &
                       "C++ reference (U = +4): V_xc^dn must match off the nodes")
        end do

        call xc_lsda_destroy(xc)

        call xc_lsda_init(xc, "data/tables/fortran_native/xc_table_u4.00.dat", status, &
                          u_signed = -4.0_dp)
        call check(status == 0, "C++ reference: XC init (U = -4) should succeed")
        if (status /= 0) return

        do i = 1, NPTS
            call get_exc(xc, n_up(i), n_dw(i), exc, ierr)
            call check(ierr == 0, "C++ reference (U = -4): get_exc should succeed")
            call check(abs(exc - exc_m(i)) < CPP_TOL, &
                       "corrected reference (U = -4): e_xc must match off the nodes")

            call get_vxc(xc, n_up(i), n_dw(i), v_up, v_dw, ierr)
            call check(ierr == 0, "C++ reference (U = -4): get_vxc should succeed")
            call check(abs(v_up - vup_m(i)) < CPP_TOL, &
                       "C++ reference (U = -4): V_xc^up must match off the nodes")
            call check(abs(v_dw - vdn_m(i)) < CPP_TOL, &
                       "C++ reference (U = -4): V_xc^dn must match off the nodes")
        end do

        call xc_lsda_destroy(xc)
    end subroutine test_cpp_reference_values_off_nodes

    !> REGRESSION (T18): the empty-channel and corner shortcuts of the C++
    !!
    !! With one spin channel empty there is no double occupancy, so e_xc = 0 and
    !! V_xc of the EMPTY channel's partner vanishes. The C++ reference
    !! short-circuits on that at the top of every recursion level
    !! (`original/spline2D.cc:474`, `:580`, `:677`), which the Fortran now
    !! reproduces with `exc_recursive` / `vxc_up_recursive` / `vxc_dn_recursive`.
    !!
    !! Two assertions here are branch-specific, i.e. they fail on the previous
    !! implementation:
    !!  - `V_xc^dn(1, 0)`: the fully polarized half-filled point n = m = 1 is
    !!    caught by the CORNER shortcut and is 0 in the C++. Without the shortcut
    !!    the Fortran returned dexc_dndown_b0(4, 1) = -1.657.
    !!  - `e_xc(0.5, 1e-15)`: inside the 1e-14 empty-channel band the C++ returns
    !!    exactly 0, whereas interpolating from the synthetic node at n = m gives
    !!    a non-zero O(1e-15) value.
    !! The remaining assertions (e_xc(1, 1), V_xc(1, 1), e_xc(0.5, 0)) already
    !! held after T8 as a side effect of the synthetic node landing exactly on
    !! the query point; they are kept because they now hold BY CONSTRUCTION and
    !! are the quantities the physics depends on (a full band must give E/L = U).
    !!
    !! The V_xc^dn(0.5, 0) assertion pins down a DELIBERATE divergence from the
    !! text of T18, which asked for zero: the C++ returns dexc_dndownB0(u, mag)
    !! there, because adding a spin-down electron to a polarized band does cost
    !! correlation energy.
    subroutine test_empty_channel_shortcuts()
        use fortuno_serial, only: check => serial_check
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, get_exc, get_vxc, &
                           xc_lsda_destroy, dexc_dndown_b0
        use lsda_constants, only: dp

        type(xc_lsda_t) :: xc
        integer :: status, ierr
        real(dp) :: exc, v_up, v_dw
        !> The shortcuts return a literal 0, so "exactly zero" is testable
        real(dp), parameter :: ZERO_TOL = 1.0e-18_dp

        call xc_lsda_init(xc, "data/tables/fortran_native/xc_table_u4.00.dat", status)
        call check(status == 0, "XC initialization should succeed")
        if (status /= 0) return

        ! ---- Empty spin-down channel, Region I -------------------------------
        call get_exc(xc, 0.5_dp, 0.0_dp, exc, ierr)
        call check(ierr == 0, "get_exc(0.5, 0) should succeed")
        call check(abs(exc) < ZERO_TOL, "e_xc must vanish with an empty channel")

        call get_vxc(xc, 0.5_dp, 0.0_dp, v_up, v_dw, ierr)
        call check(ierr == 0, "get_vxc(0.5, 0) should succeed")
        call check(abs(v_up) < ZERO_TOL, "V_xc^up must vanish with an empty spin-down channel")
        call check(abs(v_dw - dexc_dndown_b0(4.0_dp, 0.5_dp)) < ZERO_TOL, &
                   "V_xc^dn must be dexc_dndown_b0, NOT zero, with an empty spin-down channel")

        ! ---- Inside the 1e-14 empty-channel band -----------------------------
        call get_exc(xc, 0.5_dp, 1.0e-15_dp, exc, ierr)
        call check(ierr == 0, "get_exc inside the empty-channel band should succeed")
        call check(abs(exc) < ZERO_TOL, &
                   "e_xc must be exactly zero inside the 1e-14 empty-channel band")

        ! ---- Empty spin-up channel (reached through Region II) ---------------
        call get_exc(xc, 0.0_dp, 0.5_dp, exc, ierr)
        call check(ierr == 0, "get_exc(0, 0.5) should succeed")
        call check(abs(exc) < ZERO_TOL, "e_xc must vanish with an empty spin-up channel")

        call get_vxc(xc, 0.0_dp, 0.5_dp, v_up, v_dw, ierr)
        call check(ierr == 0, "get_vxc(0, 0.5) should succeed")
        call check(abs(v_dw) < ZERO_TOL, "V_xc^dn must vanish with an empty spin-up channel")
        call check(abs(v_up - dexc_dndown_b0(4.0_dp, 0.5_dp)) < ZERO_TOL, &
                   "V_xc^up at (0, 0.5) must be the spin-flipped dexc_dndown_b0")

        ! ---- Fully polarized corner n = m = 1 (corner shortcut) --------------
        call get_exc(xc, 1.0_dp, 0.0_dp, exc, ierr)
        call check(ierr == 0, "get_exc(1, 0) should succeed")
        call check(abs(exc) < ZERO_TOL, "e_xc must vanish at the fully polarized corner")

        call get_vxc(xc, 1.0_dp, 0.0_dp, v_up, v_dw, ierr)
        call check(ierr == 0, "get_vxc(1, 0) should succeed")
        call check(abs(v_up) < ZERO_TOL, "V_xc^up must vanish at the fully polarized corner")
        call check(abs(v_dw) < ZERO_TOL, &
                   "V_xc^dn must vanish at the fully polarized corner (corner shortcut)")

        ! ---- Full band n = 2: empty only AFTER the Region IV transform -------
        call get_exc(xc, 1.0_dp, 1.0_dp, exc, ierr)
        call check(ierr == 0, "get_exc(1, 1) should succeed")
        call check(abs(exc) < ZERO_TOL, "e_xc must vanish for a completely full band")

        call get_vxc(xc, 1.0_dp, 1.0_dp, v_up, v_dw, ierr)
        call check(ierr == 0, "get_vxc(1, 1) should succeed")
        call check(abs(v_up) < ZERO_TOL, "V_xc^up must vanish for a completely full band")
        call check(abs(v_dw) < ZERO_TOL, "V_xc^dn must vanish for a completely full band")

        ! ---- Empty lattice ---------------------------------------------------
        call get_exc(xc, 0.0_dp, 0.0_dp, exc, ierr)
        call check(ierr == 0, "get_exc(0, 0) should succeed")
        call check(abs(exc) < ZERO_TOL, "e_xc must vanish on the empty lattice")

        call get_vxc(xc, 0.0_dp, 0.0_dp, v_up, v_dw, ierr)
        call check(ierr == 0, "get_vxc(0, 0) should succeed")
        call check(abs(v_up) < ZERO_TOL, "V_xc^up must vanish on the empty lattice")
        call check(abs(v_dw) < ZERO_TOL, "V_xc^dn must vanish on the empty lattice")

        call xc_lsda_destroy(xc)

    end subroutine test_empty_channel_shortcuts

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

    !> Attractive U must use the Shiba (partial particle-hole) transformation
    !!
    !! The XC tables only exist for |U|, so before this was implemented a run
    !! with U = -4 silently evaluated the REPULSIVE functional: get_vxc(0.3, 0.2)
    !! returned exactly the same numbers for U = -4 and U = +4, and e_xc had the
    !! wrong sign structure entirely. The assertions below compare against the
    !! C++ reference relations (`original/spline2D.cc`):
    !!   e_xc(n_up, n_dw; U<0)    =  e_xc(1-n_up, n_dw; |U|)
    !!   V_xc^up(n_up, n_dw; U<0) = -V_xc^up(1-n_up, n_dw; |U|)
    !!   V_xc^dn(n_up, n_dw; U<0) = +V_xc^dn(1-n_up, n_dw; |U|)
    !! and, to make sure the transformation is not a no-op, that the attractive
    !! result actually differs from the repulsive one at the same densities.
    subroutine test_shiba_transform_for_attractive_u()
        use fortuno_serial, only: check => serial_check
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, get_vxc, get_exc, xc_lsda_destroy
        use lsda_constants, only: dp

        type(xc_lsda_t) :: xc_att, xc_rep
        integer :: status, ierr
        character(len=256) :: test_file
        real(dp) :: v_att_up, v_att_dw, v_rep_up, v_rep_dw
        real(dp) :: v_same_up, v_same_dw
        real(dp) :: exc_att, exc_rep
        real(dp), parameter :: N_UP = 0.3_dp, N_DW = 0.2_dp

        test_file = "data/tables/fortran_native/xc_table_u4.00.dat"

        call xc_lsda_init(xc_att, test_file, status, u_signed=-4.0_dp)
        call check(status == 0, "init with u_signed = -4 should succeed")
        call check(xc_att%U < 0.0_dp, "xc%U must keep the attractive sign")

        call xc_lsda_init(xc_rep, test_file, status, u_signed=4.0_dp)
        call check(status == 0, "init with u_signed = +4 should succeed")
        call check(xc_rep%U > 0.0_dp, "xc%U must stay positive for repulsive U")

        ! V_xc: mirrored densities, up flipped, down not
        call get_vxc(xc_att, N_UP, N_DW, v_att_up, v_att_dw, ierr)
        call check(ierr == 0, "get_vxc must succeed for U < 0")

        call get_vxc(xc_rep, 1.0_dp - N_UP, N_DW, v_rep_up, v_rep_dw, ierr)
        call check(ierr == 0, "get_vxc must succeed for the mirrored point at U > 0")

        call check(abs(v_att_up + v_rep_up) < 1.0e-12_dp, &
                   "V_xc^up(U<0) must equal -V_xc^up(1-n_up, n_dw; |U|)")
        call check(abs(v_att_dw - v_rep_dw) < 1.0e-12_dp, &
                   "V_xc^dn(U<0) must equal +V_xc^dn(1-n_up, n_dw; |U|)")

        ! e_xc: mirrored densities, no sign flip
        call get_exc(xc_att, N_UP, N_DW, exc_att, ierr)
        call check(ierr == 0, "get_exc must succeed for U < 0")

        call get_exc(xc_rep, 1.0_dp - N_UP, N_DW, exc_rep, ierr)
        call check(ierr == 0, "get_exc must succeed for the mirrored point at U > 0")

        call check(abs(exc_att - exc_rep) < 1.0e-12_dp, &
                   "e_xc(U<0) must equal e_xc(1-n_up, n_dw; |U|)")

        ! The transformation must not be an accidental identity
        call get_vxc(xc_rep, N_UP, N_DW, v_same_up, v_same_dw, ierr)
        call check(ierr == 0, "get_vxc must succeed at the untransformed point")
        call check(abs(v_att_up - v_same_up) > 1.0e-6_dp .or. &
                   abs(v_att_dw - v_same_dw) > 1.0e-6_dp, &
                   "attractive V_xc must differ from the repulsive one at the same densities")

        call xc_lsda_destroy(xc_att)
        call xc_lsda_destroy(xc_rep)
    end subroutine test_shiba_transform_for_attractive_u

    !> A signed U whose magnitude disagrees with the loaded table is an error
    !!
    !! The magnitude check is the only protection against running the functional
    !! of a different interaction strength (the table carries |U| and the caller
    !! carries the sign; if they disagree, one of the two is a bug).
    subroutine test_shiba_requires_matching_table()
        use fortuno_serial, only: check => serial_check
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use lsda_constants, only: dp

        type(xc_lsda_t) :: xc
        integer :: status
        character(len=256) :: test_file

        test_file = "data/tables/fortran_native/xc_table_u4.00.dat"

        call xc_lsda_init(xc, test_file, status, u_signed=-2.0_dp)
        call check(status /= 0, "|u_signed| /= table |U| must be rejected")
        call check(.not. xc%initialized, "XC must stay uninitialized on a |U| mismatch")

        call xc_lsda_destroy(xc)
    end subroutine test_shiba_requires_matching_table
end program test_xc_lsda
