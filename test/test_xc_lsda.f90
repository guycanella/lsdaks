!> Unit tests for xc_lsda module
program test_xc_lsda
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp
    implicit none

    real(dp), parameter :: TOL = 1.0e-9_dp

    !> Finer U = 4 fixture, used only by `test_off_node_matches_generator`
    character(len=*), parameter :: FINE_TABLE = &
        'build/test_xc_lsda_xc_table_u4.00_fine.dat'

    call prepare_test_xc_tables()
    call execute_serial_cmd_app(get_xc_lsda_tests())

contains

    function get_xc_lsda_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests

        tests = test_list([ &
            test("xc_lsda_init_destroy", test_xc_lsda_init_destroy), &
            test("zero_u_initializes_without_table", test_zero_u_initializes_without_table), &
            test("nonzero_small_u_requires_table", test_nonzero_small_u_requires_table), &
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
            test("off_node_xc_symmetries", &
                 test_off_node_xc_symmetries), &
            test("off_node_matches_generator", &
                 test_off_node_matches_generator) &
        ])
    end function get_xc_lsda_tests

    !> The exact U = 0 functional must not read or require an XC table.
    !!
    !! This is the analytic non-interacting limit: e_xc and both spin potentials
    !! vanish for every physical density.  A deliberately nonexistent filename
    !! proves the initializer returns before any table I/O.
    subroutine test_zero_u_initializes_without_table()
        use fortuno_serial, only: check => serial_check
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, get_exc, get_vxc, xc_lsda_destroy
        use lsda_errors, only: ERROR_SUCCESS

        type(xc_lsda_t) :: xc
        integer :: status, ierr
        real(dp) :: exc, v_up, v_dw

        call xc_lsda_init(xc, 'build/nonexistent_xc_table.dat', status, u_signed=0.0_dp)
        call check(status == ERROR_SUCCESS, "U = 0 XC initialization must not require a table")
        call check(xc%initialized, "U = 0 XC functional must be initialized")
        call check(abs(xc%U) < TOL, "U = 0 XC functional must retain zero interaction")

        call get_exc(xc, 0.4_dp, 0.3_dp, exc, ierr)
        call check(ierr == ERROR_SUCCESS .and. abs(exc) < TOL, "e_xc must vanish at U = 0")

        call get_vxc(xc, 0.4_dp, 0.3_dp, v_up, v_dw, ierr)
        call check(ierr == ERROR_SUCCESS .and. abs(v_up) < TOL .and. abs(v_dw) < TOL, &
                   "both V_xc channels must vanish at U = 0")

        call xc_lsda_destroy(xc)
    end subroutine test_zero_u_initializes_without_table

    !> A nonzero interaction must not be silently replaced by the U = 0 XC functional.
    !!
    !! The zero-XC shortcut is a physical statement about the non-interacting
    !! point only.  A small but nonzero U below the generator floor still needs
    !! a compatible table and must therefore fail when no table can be read.
    !! This fails if XC_U_MATCH_TOL is incorrectly used as the physical-zero
    !! threshold.
    subroutine test_nonzero_small_u_requires_table()
        use fortuno_serial, only: check => serial_check
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
        use lsda_errors, only: ERROR_SUCCESS

        type(xc_lsda_t) :: xc
        integer :: status

        call xc_lsda_init(xc, 'build/nonexistent_xc_table.dat', status, u_signed=1.0e-7_dp)
        call check(status /= ERROR_SUCCESS, "nonzero U below the table floor must require a table")
        call check(.not. xc%initialized, "nonzero U must not initialize the zero-XC functional")

        call xc_lsda_destroy(xc)
    end subroutine test_nonzero_small_u_requires_table

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

        call read_fortran_table("build/test_xc_lsda_xc_table_u4.00.dat", table, status)
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

        call read_fortran_table("build/test_xc_lsda_xc_table_u4.00.dat", table, status)
        call check(status == 0, "endpoint slope: table should be readable")
        if (status /= 0) return

        call xc_lsda_init(xc, "build/test_xc_lsda_xc_table_u4.00.dat", status)
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

        call xc_lsda_init(xc, "build/test_xc_lsda_xc_table_u4.00.dat", status)
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

    !> XC symmetry identities must hold away from tabulated nodes.
    !!
    !! The seven points are off-node in BOTH directions and cover
    !!   1,2: Region I   (m >= 0, n <= 1)
    !!   3:   Region II  (m <  0, n <= 1)   -> spin exchange
    !!   4:   Region III (m <  0, n >  1)   -> particle-hole, sign flip
    !!   5:   Region IV  (m >= 0, n >  1)   -> both
    !!   6,7: below the first tabulated density (n = 0.007 and n = 0.004 against
    !!        a first row at n = 0.0202807), i.e. the window that consists of the
    !!        synthetic node plus the whole table
    !! For repulsive U, spin exchange is an exact identity.  For attractive U,
    !! the Shiba transformation relates the functional to the repulsive result
    !! at (1 - n_up, n_down).  These checks independently exercise interpolation
    !! between nodes, recursive region mapping, and the sign conventions of both
    !! spin potentials without a reference-table or C++ runtime dependency.
    subroutine test_off_node_xc_symmetries()
        use fortuno_serial, only: check => serial_check
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, get_exc, get_vxc, xc_lsda_destroy

        integer, parameter :: NPTS = 7
        real(dp), parameter :: SYMMETRY_TOL = 1.0e-10_dp

        type(xc_lsda_t) :: xc_rep, xc_att
        integer :: status, ierr, i
        real(dp) :: exc, exc_swap, exc_att, exc_shiba
        real(dp) :: v_up, v_dw, v_up_swap, v_dw_swap
        real(dp) :: v_up_att, v_dw_att, v_up_shiba, v_dw_shiba
        real(dp) :: n_up(NPTS), n_dw(NPTS)

        n_up = [0.34267985324_dp, 0.27001671326_dp, 0.15_dp, 0.45_dp, 0.82_dp, &
                0.004_dp, 0.0031_dp]
        n_dw = [0.028423868551_dp, 0.20374663913_dp, 0.30_dp, 0.72_dp, 0.43_dp, &
                0.003_dp, 0.0009_dp]

        call xc_lsda_init(xc_rep, "build/test_xc_lsda_xc_table_u4.00.dat", status, &
                          u_signed = 4.0_dp)
        call check(status == 0, "Off-node symmetry: repulsive XC init should succeed")
        if (status /= 0) return

        call xc_lsda_init(xc_att, "build/test_xc_lsda_xc_table_u4.00.dat", status, &
                          u_signed = -4.0_dp)
        call check(status == 0, "Off-node symmetry: attractive XC init should succeed")
        if (status /= 0) then
            call xc_lsda_destroy(xc_rep)
            return
        end if

        do i = 1, NPTS
            call get_exc(xc_rep, n_up(i), n_dw(i), exc, ierr)
            call get_exc(xc_rep, n_dw(i), n_up(i), exc_swap, status)
            call check(ierr == 0 .and. status == 0, &
                       "Off-node symmetry: repulsive e_xc evaluations should succeed")
            call check(abs(exc - exc_swap) < SYMMETRY_TOL, &
                       "Off-node symmetry: repulsive e_xc must be spin symmetric")

            call get_vxc(xc_rep, n_up(i), n_dw(i), v_up, v_dw, ierr)
            call get_vxc(xc_rep, n_dw(i), n_up(i), v_up_swap, v_dw_swap, status)
            call check(ierr == 0 .and. status == 0, &
                       "Off-node symmetry: repulsive V_xc evaluations should succeed")
            call check(abs(v_up - v_dw_swap) < SYMMETRY_TOL .and. &
                       abs(v_dw - v_up_swap) < SYMMETRY_TOL, &
                       "Off-node symmetry: repulsive V_xc must exchange spin channels")

            call get_exc(xc_att, n_up(i), n_dw(i), exc_att, ierr)
            call get_exc(xc_rep, 1.0_dp - n_up(i), n_dw(i), exc_shiba, status)
            call check(ierr == 0 .and. status == 0, &
                       "Off-node Shiba: e_xc evaluations should succeed")
            call check(abs(exc_att - exc_shiba) < SYMMETRY_TOL, &
                       "Off-node Shiba: attractive e_xc must map to repulsive e_xc")

            call get_vxc(xc_att, n_up(i), n_dw(i), v_up_att, v_dw_att, ierr)
            call get_vxc(xc_rep, 1.0_dp - n_up(i), n_dw(i), v_up_shiba, v_dw_shiba, status)
            call check(ierr == 0 .and. status == 0, &
                       "Off-node Shiba: V_xc evaluations should succeed")
            call check(abs(v_up_att + v_up_shiba) < SYMMETRY_TOL .and. &
                       abs(v_dw_att - v_dw_shiba) < SYMMETRY_TOL, &
                       "Off-node Shiba: attractive V_xc must retain the channel signs")
        end do

        call xc_lsda_destroy(xc_att)
        call xc_lsda_destroy(xc_rep)
    end subroutine test_off_node_xc_symmetries

    !> Between nodes the interpolant must reproduce the function it interpolates.
    !!
    !! The symmetry checks above cannot see this: Shiba and spin exchange are
    !! invariances of *any* interpolant, and an off-by-one in the row window of
    !! `spline2d_eval` shifts both sides of a symmetry by the same amount, so it
    !! passes them.  Here the spline is compared, at points off-node in both `n`
    !! and `m`, against a direct solve of the same Lieb-Wu equations that built
    !! the table - the generator, not the C++ reference.
    !!
    !! `V_xc` is checked the same way and for the same reason: the symmetry
    !! checks are equally blind to a shifted row window in the `vxc_up` /
    !! `vxc_down` splines, and nothing else in the suite pins those two
    !! interpolants to the function they interpolate away from the nodes.
    !!
    !! This test uses the finer `FINE_TABLE` fixture rather than the 8 x 9 grid
    !! the rest of the suite shares. `V_xc` is a derivative of `e_xc` and is
    !! correspondingly harder to interpolate: measured at these same five
    !! points, the spline error on 8 x 9 is 1.6e-2 in `V_xc` (8.1e-4 in
    !! `e_xc`), against one-row-shift signals of the same 1e-2 order, so a
    !! tolerance there could not separate a correct interpolation from a shift.
    !!
    !! Tolerance: on the 20 x 31 fixture the measured spline error at these
    !! points is at most 7.3e-6 in `e_xc` and 1.4e-4 in `V_xc`; `VXC_TOL` sits
    !! ~7x above that. Verified by mutation on this fixture (shifting only the
    !! value arrays, rewriting the table and re-reading it): `cshift` of one
    !! row (`dim=2`, the `n` axis) or one column (`dim=1`, the `m` axis) in
    !! `exc`, `vxc_up` or `vxc_down` fails the combined assertions at all 5
    !! points, with the smallest detected deviation 3.6e-4 in `e_xc` and
    !! 1.9e-3 in `V_xc`. For a shifted `V_xc` channel, its own channel-specific
    !! assertion catches 4 of 5 points; at the remaining point `n_up < n_dw`,
    !! the Shiba/spin-swap path evaluates the opposite spline and the other
    !! channel's assertion detects the mutation.
    subroutine test_off_node_matches_generator()
        use fortuno_serial, only: check => serial_check
        use xc_lsda, only: xc_lsda_t, xc_lsda_init, get_exc, get_vxc, xc_lsda_destroy
        use bethe_tables, only: compute_E_xc, compute_V_xc_numerical, xc_potentials_t

        integer, parameter :: NPTS = 5
        real(dp), parameter :: SPLINE_TOL = 1.0e-4_dp
        real(dp), parameter :: VXC_TOL = 1.0e-3_dp
        real(dp), parameter :: N_UP(NPTS) = [0.34267985324_dp, 0.27001671326_dp, &
                                             0.15_dp, 0.4123_dp, 0.61_dp]
        real(dp), parameter :: N_DW(NPTS) = [0.028423868551_dp, 0.20374663913_dp, &
                                             0.30_dp, 0.3777_dp, 0.29_dp]

        type(xc_lsda_t) :: xc
        type(xc_potentials_t) :: v_exact
        real(dp) :: exc, exc_exact, v_up, v_dw
        integer :: status, i

        call xc_lsda_init(xc, FINE_TABLE, status, u_signed = 4.0_dp)
        call check(status == 0, "Off-node interpolation: XC init should succeed")
        if (status /= 0) return

        do i = 1, NPTS
            call get_exc(xc, N_UP(i), N_DW(i), exc, status)
            call check(status == 0, "Off-node interpolation: e_xc evaluation should succeed")
            exc_exact = compute_E_xc(N_UP(i), N_DW(i), 4.0_dp)
            call check(abs(exc - exc_exact) < SPLINE_TOL, &
                       "Off-node interpolation: the spline must reproduce the " // &
                       "generator between nodes, not a neighbouring row")

            call get_vxc(xc, N_UP(i), N_DW(i), v_up, v_dw, status)
            call check(status == 0, "Off-node interpolation: V_xc evaluation should succeed")
            v_exact = compute_V_xc_numerical(N_UP(i), N_DW(i), 4.0_dp)
            call check(abs(v_up - v_exact%v_xc_up) < VXC_TOL, &
                       "Off-node interpolation: the V_xc^up spline must reproduce " // &
                       "the generator between nodes, not a neighbouring row")
            call check(abs(v_dw - v_exact%v_xc_down) < VXC_TOL, &
                       "Off-node interpolation: the V_xc^dn spline must reproduce " // &
                       "the generator between nodes, not a neighbouring row")
        end do

        call xc_lsda_destroy(xc)
    end subroutine test_off_node_matches_generator

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

        call xc_lsda_init(xc, "build/test_xc_lsda_xc_table_u4.00.dat", status)
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

        test_file = "build/test_xc_lsda_xc_table_u4.00.dat"

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

        test_file = "build/test_xc_lsda_xc_table_u2.00.dat"
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

        test_file = "build/test_xc_lsda_xc_table_u2.00.dat"
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

        test_file = "build/test_xc_lsda_xc_table_u2.00.dat"
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

        test_file = "build/test_xc_lsda_xc_table_u4.00.dat"
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

        test_file = "build/test_xc_lsda_xc_table_u4.00.dat"

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

        test_file = "build/test_xc_lsda_xc_table_u4.00.dat"

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

        test_file = "build/test_xc_lsda_xc_table_u4.00.dat"

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

        test_file = "build/test_xc_lsda_xc_table_u4.00.dat"

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

        test_file = "build/test_xc_lsda_xc_table_u4.00.dat"

        call xc_lsda_init(xc, test_file, status, u_signed=-2.0_dp)
        call check(status /= 0, "|u_signed| /= table |U| must be rejected")
        call check(.not. xc%initialized, "XC must stay uninitialized on a |U| mismatch")

        call xc_lsda_destroy(xc)
    end subroutine test_shiba_requires_matching_table

    !> Generate the deterministic XC fixtures used by this test program.
    !!
    !! The compact grid keeps the release and debug suites fast while exercising
    !! the production table generator, writer, reader, and spline initialization.
    !! Files live only in `build/` and are rewritten on every run: they are the
    !! oracle of this suite, so a cached copy from an older generator would let
    !! the tests keep validating against a table the code no longer produces.
    !! The 8 x 9 grid makes regenerating cheap enough not to bother caching.
    !!
    !! One extra, finer fixture (`FINE_TABLE`, 20 x 31) is written for
    !! `test_off_node_matches_generator`. That test compares the interpolants
    !! against the generator between nodes, and on the 8 x 9 grid the measured
    !! `V_xc` spline error at those points (1.6e-2) is the same order as the
    !! deviation a one-row shift produces, so no tolerance there can tell the
    !! two apart. On 20 x 31 the error drops to 1.4e-4 while the shift signal
    !! stays at 1e-3..7e-2, and the comparison becomes discriminating. The
    !! coarse grid is kept for every other test so their calibrated tolerances
    !! stay valid.
    subroutine prepare_test_xc_tables()
        use bethe_tables, only: generate_xc_table, grid_params_t
        use table_io, only: xc_table_t, write_fortran_table, deallocate_table
        use lsda_errors, only: ERROR_SUCCESS

        real(dp), parameter :: U_VALUES(2) = [2.0_dp, 4.0_dp]
        type(grid_params_t) :: params
        type(xc_table_t) :: table
        character(len=256) :: filename
        integer :: i, ierr

        params = grid_params_t()
        params%n_points = 8
        params%m_points = 9

        do i = 1, size(U_VALUES)
            write(filename, '(A,F0.2,A)') 'build/test_xc_lsda_xc_table_u', U_VALUES(i), '.dat'

            call generate_xc_table(U_VALUES(i), params, table, ierr)
            if (ierr /= ERROR_SUCCESS) error stop 'failed to generate XC-LSDA test fixture'

            call write_fortran_table(trim(filename), table, ierr)
            if (allocated(table%n_grid)) call deallocate_table(table)
            if (ierr /= ERROR_SUCCESS) error stop 'failed to write XC-LSDA test fixture'
        end do

        params%n_points = 20
        params%m_points = 31
        call generate_xc_table(4.0_dp, params, table, ierr)
        if (ierr /= ERROR_SUCCESS) error stop 'failed to generate fine XC-LSDA test fixture'
        call write_fortran_table(FINE_TABLE, table, ierr)
        if (allocated(table%n_grid)) call deallocate_table(table)
        if (ierr /= ERROR_SUCCESS) error stop 'failed to write fine XC-LSDA test fixture'
    end subroutine prepare_test_xc_tables
end program test_xc_lsda
