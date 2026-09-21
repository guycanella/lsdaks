!> Unit tests for bethe_tables module
program test_bethe_tables
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
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
            test("compute_V_xc_polarized_corner", test_compute_V_xc_polarized_corner), &
            test("attractive_shiba_internal_consistency_u4", test_attractive_shiba_u4), &
            test("vxc_edge_identity_only_on_edge_u4", test_vxc_edge_identity_only_on_edge), &
            test("graded_grids", test_graded_grids), &
            test("table_generator_u_floor", test_table_generator_u_floor), &
            test("half_filling_and_edge_guards_u4", test_half_filling_and_edge_guards), &
            test("compute_E_xc_above_half_filling_u4", test_compute_E_xc_above_half_filling_u4), &
            test("compute_E_xc_particle_hole_u4", test_compute_E_xc_particle_hole_u4), &
            test("compute_E_xc_spin_exchange_u4", test_compute_E_xc_spin_exchange_u4), &
            test("compute_V_xc_symmetric", test_compute_V_xc_symmetric), &
            test("grid_params_defaults", test_grid_params_defaults), &
            test("invalid_density_grid_bounds", test_invalid_density_grid_bounds), &
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

        n_up = 0.3_dp
        n_dw = 0.2_dp
        U = 0.0_dp
        
        E_xc = compute_E_xc(n_up, n_dw, U)
        print *, "E_xc for U=0:", E_xc
        
        ! For U=0, E_xc should be zero (E_BA = E0)
        call check(abs(E_xc) < 1.0e-6_dp, "E_xc should be ~0 for U=0")
        
    end subroutine test_compute_E_xc_U0

    !> Regression test for the U=4 half-filled thermodynamic-limit solution.
    !!
    !! The reference is the C++ U=4 table value after the Hartree contribution
    !! has been removed.  The reference file prints this row with ~9 significant
    !! digits, so 1e-8 is the tightest meaningful bound against it.
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

        E_xc = compute_E_xc(0.5_dp, 0.5_dp, 4.0_dp)
        call check(abs(E_xc - CPP_EXC) < 1.0e-8_dp, &
                   "U=4 half-filled E_xc must match the C++ reference table")

        params%n_min = 1.0_dp
        params%n_max = 1.0_dp
        params%n_points = 1
        params%m_points = 1
        call generate_xc_table(4.0_dp, params, table, status)
        call check(status == 0, "U=4 reference table point should generate")
        call check(abs(table%exc(1, 1) - CPP_EXC) < 1.0e-8_dp, &
                   "Generated U=4 table must match the C++ reference table")

        deallocate(table%n_grid, table%m_grid, table%exc, table%vxc_up, table%vxc_down)
    end subroutine test_compute_E_xc_half_filling_u4

    !> Regression test for the U=4 potentials at half filling.
    !!
    !! `e_xc` has a cusp at `n = 1` (Mott gap) and is symmetric about
    !! `n_up = 0.5` at fixed `n_dn = 0.5` by particle-hole symmetry, so a
    !! central difference across half filling returns exactly zero.  The
    !! tabulated value is the one-sided limit from `n < 1`, which is what this
    !! test pins: it fails with 0 instead of -0.6434 if the stencil is allowed
    !! to step above half filling.
    subroutine test_compute_V_xc_half_filling_u4()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: compute_V_xc_numerical, xc_potentials_t
        use lsda_constants, only: dp

        real(dp), parameter :: CPP_VXC = -0.643363524_dp
        type(xc_potentials_t) :: v_xc

        v_xc = compute_V_xc_numerical(0.5_dp, 0.5_dp, 4.0_dp)

        call check(ieee_is_finite(v_xc%v_xc_up), &
                   "U=4 half-filled V_xc_up must be finite")
        call check(ieee_is_finite(v_xc%v_xc_down), &
                   "U=4 half-filled V_xc_down must be finite")
        call check(abs(v_xc%v_xc_up - v_xc%v_xc_down) < 1.0e-14_dp, &
                   "U=4 half-filled V_xc must preserve spin symmetry")
        call check(abs(v_xc%v_xc_up - CPP_VXC) < 1.0e-6_dp, &
                   "U=4 half-filled V_xc must be the one-sided limit from n < 1")
    end subroutine test_compute_V_xc_half_filling_u4

    !> Analytic regression test for the fully polarized corner at several U.
    !!
    !! At `n = 1`, `m = 1` the only stencils that stay inside the physical
    !! triangle are the two edges: `e_xc == 0` along `m = n` forces
    !! `V_up = 0`, and the `m` derivative along `n = 1` gives
    !! `V_dn = 4 - sqrt(U^2 + 16)` for every repulsive U.
    subroutine test_compute_V_xc_polarized_corner()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: compute_E_xc, compute_V_xc_numerical, xc_potentials_t
        use lsda_constants, only: dp

        !> Closed form of the corner value. The previous U = 4-only instance,
        !! `-4 (sqrt(2) - 1)`, coincides with several plausible expressions;
        !! sweeping U makes the analytic dependence independently auditable.
        !!
        !! The tolerance cannot go all the way down to the working precision:
        !! `V_dn` here is a finite-difference derivative along the `n = 1` edge,
        !! and the measured deviation of the stencil from the closed form is
        !! 2.7e-8.  1e-7 is therefore the tightest honest bound, still an order
        !! of magnitude below the 1e-6 this test used to accept.
        real(dp), parameter :: U_VALUES(4) = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp]
        type(xc_potentials_t) :: v_xc
        real(dp) :: e_xc, vxc_dn_exact
        integer :: i

        do i = 1, size(U_VALUES)
            e_xc = compute_E_xc(1.0_dp, 0.0_dp, U_VALUES(i))
            v_xc = compute_V_xc_numerical(1.0_dp, 0.0_dp, U_VALUES(i))
            vxc_dn_exact = 4.0_dp - sqrt(U_VALUES(i)**2 + 16.0_dp)

            call check(abs(e_xc) < 1.0e-14_dp, "e_xc must vanish at the polarized corner")
            call check(abs(v_xc%v_xc_up) < 1.0e-14_dp, &
                       "V_xc_up must vanish at the polarized corner")
            call check(abs(v_xc%v_xc_down - vxc_dn_exact) < 1.0e-7_dp, &
                       "V_xc_down at the polarized corner must match its analytic form")
        end do
    end subroutine test_compute_V_xc_polarized_corner

    !> **Internal consistency only** of the attractive Shiba branch.
    !!
    !! `e_xc(n_up, n_dn; -U) = e_xc(1 - n_up, n_dn; U)`, with `V_xc_up`
    !! changing sign and `V_xc_down` keeping it.  Both sides of this test go
    !! through the same `shiba_map`, so it pins the convention the module uses
    !! internally and nothing more: flipping the channel that changes sign in
    !! **both** the generator and `xc_lsda` would leave this test passing.
    subroutine test_attractive_shiba_u4()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: compute_E_xc, compute_V_xc_numerical, xc_potentials_t
        use lsda_constants, only: dp

        type(xc_potentials_t) :: v_att, v_rep
        real(dp) :: e_att, e_rep

        e_att = compute_E_xc(0.3_dp, 0.2_dp, -4.0_dp)
        e_rep = compute_E_xc(0.7_dp, 0.2_dp, 4.0_dp)
        call check(abs(e_att - e_rep) < 1.0e-14_dp, &
                   "attractive e_xc must equal the Shiba-mapped repulsive one")

        v_att = compute_V_xc_numerical(0.3_dp, 0.2_dp, -4.0_dp)
        v_rep = compute_V_xc_numerical(0.7_dp, 0.2_dp, 4.0_dp)
        call check(abs(v_att%v_xc_up + v_rep%v_xc_up) < 1.0e-14_dp, &
                   "attractive V_xc_up must be minus the Shiba-mapped one")
        call check(abs(v_att%v_xc_down - v_rep%v_xc_down) < 1.0e-14_dp, &
                   "attractive V_xc_down must keep the Shiba-mapped sign")
    end subroutine test_attractive_shiba_u4

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

        e_xc = compute_E_xc(0.75_dp, 0.75_dp, 4.0_dp)
        v_xc = compute_V_xc_numerical(0.75_dp, 0.75_dp, 4.0_dp)

        call check(ieee_is_finite(e_xc), "Above-half-filled U=4 E_xc must be finite")
        call check(ieee_is_finite(v_xc%v_xc_up), &
                   "Above-half-filled U=4 V_xc_up must be finite")
        call check(ieee_is_finite(v_xc%v_xc_down), &
                   "Above-half-filled U=4 V_xc_down must be finite")
    end subroutine test_compute_E_xc_above_half_filling_u4

    !> Regression test for particle-hole symmetry of the XC energy.
    !!
    !! The Lieb-Wu solve above half filling is performed at complementary
    !! densities, so the Hartree term must use those same mapped densities.
    subroutine test_compute_E_xc_particle_hole_u4()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: compute_E_xc
        use lsda_constants, only: dp

        real(dp) :: e_xc, e_xc_complement

        e_xc = compute_E_xc(0.75_dp, 0.75_dp, 4.0_dp)
        e_xc_complement = compute_E_xc(0.25_dp, 0.25_dp, 4.0_dp)

        call check(abs(e_xc - e_xc_complement) < 1.0e-12_dp, &
                   "E_xc must obey particle-hole symmetry above half filling")
    end subroutine test_compute_E_xc_particle_hole_u4

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

        e_xc_up = compute_E_xc(0.1_dp, 0.7_dp, 4.0_dp)
        e_xc_down = compute_E_xc(0.7_dp, 0.1_dp, 4.0_dp)

        call check(ieee_is_finite(e_xc_up), &
                   "Negative-m state must not enter an invalid Bethe sector")
        call check(ieee_is_finite(e_xc_down), &
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
        type(xc_potentials_t) :: v_xc

        n_up = 0.25_dp
        n_dw = 0.25_dp  ! Symmetric
        U = 1.0_dp

        v_xc = compute_V_xc_numerical(n_up, n_dw, U)

        ! For symmetric case, potentials should be equal
        call check(abs(v_xc%v_xc_up - v_xc%v_xc_down) < 1.0e-6_dp, &
                   "V_xc_up should equal V_xc_down for symmetric case")
        
    end subroutine test_compute_V_xc_symmetric

    !> Test grid_params_t default values
    !!
    !! The grid sizes and grading exponents are not free style choices: they
    !! are the operating point measured in phase 4.5 as the post-spline
    !! deviation from the reference U=4 table (2.7e-7 in e_xc, 7.7e-5 in V_xc).
    !! Changing them without repeating that measurement silently degrades every
    !! table the generator produces, so they are pinned here.
    subroutine test_grid_params_defaults()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: grid_params_t
        use lsda_constants, only: dp
        
        type(grid_params_t) :: params
        
        ! Check defaults
        call check(ABS(params%n_min - 0.02_dp) < TOL, "Default n_min should be 0.02")
        call check(ABS(params%n_max - 1.0_dp) < TOL, "Default n_max should be 1.0")
        call check(params%n_points == 75, "Default n_points should be 75")
        call check(params%m_points == 202, "Default m_points should be 202")
        call check(params%m_grade > 0.0_dp .and. params%m_grade < 1.0_dp, &
                   "Default m axis must be graded, not uniform")
        call check(params%n_grade_low > 1.0_dp .and. params%n_grade_high > 1.0_dp, &
                   "Default n axis must be compressed at both ends")
        call check(params%m_frac_min > 0.0_dp .and. params%m_frac_min < 1.0e-3_dp, &
                   "Default m axis must reach at least three decades below m = n")
        call check(params%quad%n_k >= 16, "Default k quadrature order must be usable")
        call check(params%quad%n_lambda >= 8, "Default Lambda quadrature order must be usable")
        call check(params%quad%tol > 0.0_dp, "Default inversion tolerance must be positive")
        call check(params%delta_n > 0.0_dp, "Default finite-difference step must be positive")
        
    end subroutine test_grid_params_defaults

    !> The table generator must reject nonphysical density axes before solving.
    subroutine test_invalid_density_grid_bounds()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: grid_params_t, generate_xc_table
        use table_io, only: xc_table_t
        use lsda_errors, only: ERROR_INVALID_INPUT

        type(grid_params_t) :: params
        type(xc_table_t) :: table
        integer :: status

        params%n_min = 0.0_dp
        call generate_xc_table(4.0_dp, params, table, status)
        call check(status == ERROR_INVALID_INPUT, "n_min = 0 must be rejected")

        params = grid_params_t()
        params%n_min = 0.8_dp
        params%n_max = 0.7_dp
        call generate_xc_table(4.0_dp, params, table, status)
        call check(status == ERROR_INVALID_INPUT, "a reversed density axis must be rejected")

        params = grid_params_t()
        params%n_max = 1.0_dp + 1.0e-9_dp
        call generate_xc_table(4.0_dp, params, table, status)
        call check(status == ERROR_INVALID_INPUT, &
                   "a table density axis above half filling must be rejected")

        params = grid_params_t()
        params%delta_n = 0.0_dp
        call generate_xc_table(4.0_dp, params, table, status)
        call check(status == ERROR_INVALID_INPUT, "delta_n = 0 must be rejected")

        params = grid_params_t()
        params%m_grade = ieee_value(0.0_dp, ieee_quiet_nan)
        call generate_xc_table(4.0_dp, params, table, status)
        call check(status == ERROR_INVALID_INPUT, "NaN m_grade must be rejected")

        params = grid_params_t()
        params%quad%tol = ieee_value(0.0_dp, ieee_quiet_nan)
        call generate_xc_table(4.0_dp, params, table, status)
        call check(status == ERROR_INVALID_INPUT, "NaN quadrature tolerance must be rejected")

        params = grid_params_t()
        call generate_xc_table(ieee_value(0.0_dp, ieee_quiet_nan), params, table, status)
        call check(status == ERROR_INVALID_INPUT, "NaN U must be rejected")
    end subroutine test_invalid_density_grid_bounds

    !> Regression test: the `de/dn = -de/dm` identity holds only on `m = n`.
    !!
    !! `e_xc` vanishes identically along the fully polarized edge, so there the
    !! directional derivative along the edge is zero and `V_xc_up = 0`.  The
    !! identity used to be selected by the finite-difference step
    !! (`m + 2h > n`) instead of by the distance from the edge, so a large
    !! `delta_n` made **interior** points inherit the edge identity: they came
    !! out with `V_xc_up` exactly zero and an O(n - m) error in `V_xc_down`,
    !! with no warning.  Here `n = 0.5`, `m = 0.46`, `delta_n = 0.05` gives
    !! `h = 0.025` and `m + 2h = 0.51 > n`, which is precisely the case that
    !! used to be misclassified.
    subroutine test_vxc_edge_identity_only_on_edge()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: compute_V_xc_numerical, xc_potentials_t
        use lsda_constants, only: dp

        type(xc_potentials_t) :: v_big_h, v_small_h, v_edge

        ! n = 0.5, m = 0.46 -> n_up = 0.48, n_dn = 0.02
        v_big_h = compute_V_xc_numerical(0.48_dp, 0.02_dp, 4.0_dp, delta_n=0.05_dp)
        v_small_h = compute_V_xc_numerical(0.48_dp, 0.02_dp, 4.0_dp, delta_n=1.0e-4_dp)

        call check(abs(v_big_h%v_xc_up) > 1.0e-3_dp, &
                   "an interior point must not inherit V_xc_up = 0 from the m = n edge")
        call check(abs(v_big_h%v_xc_up - v_small_h%v_xc_up) < 5.0e-3_dp, &
                   "the interior V_xc_up must not depend on the stencil step")
        call check(abs(v_big_h%v_xc_down - v_small_h%v_xc_down) < 5.0e-3_dp, &
                   "the interior V_xc_down must not depend on the stencil step")

        ! On the edge itself the identity must still be used.
        v_edge = compute_V_xc_numerical(0.5_dp, 0.0_dp, 4.0_dp, delta_n=0.05_dp)
        call check(abs(v_edge%v_xc_up) < 1.0e-14_dp, &
                   "V_xc_up must vanish exactly on the fully polarized edge")
    end subroutine test_vxc_edge_identity_only_on_edge

    !> The generated axes must be graded, monotone and hit their end points.
    subroutine test_graded_grids()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: magnetization_grid, density_grid
        use lsda_constants, only: dp

        real(dp) :: m(50), n(40)
        integer :: k
        logical :: monotone

        call magnetization_grid(0.8_dp, 1.0e-5_dp, 0.4_dp, m)
        call check(abs(m(1)) < 1.0e-300_dp, "m axis must start at 0")
        call check(abs(m(2) - 0.8_dp * 1.0e-5_dp) < 1.0e-14_dp, &
                   "m axis must start its ladder at m_frac_min * n")
        call check(abs(m(50) - 0.8_dp) < 1.0e-14_dp, "m axis must end at m = n")
        monotone = .true.
        do k = 2, 50
            if (m(k) <= m(k - 1)) monotone = .false.
        end do
        call check(monotone, "m axis must be strictly increasing")
        call check(m(3) - m(2) < 0.05_dp * (m(50) - m(49)), &
                   "m axis must be far finer next to m = 0 than next to m = n")

        call density_grid(0.02_dp, 1.0_dp, 1.1_dp, 1.45_dp, n)
        call check(abs(n(1) - 0.02_dp) < 1.0e-14_dp, "n axis must start at n_min")
        call check(abs(n(40) - 1.0_dp) < 1.0e-14_dp, "n axis must end at n_max")
        monotone = .true.
        do k = 2, 40
            if (n(k) <= n(k - 1)) monotone = .false.
        end do
        call check(monotone, "n axis must be strictly increasing")
        call check(n(40) - n(39) < 0.5_dp * (n(21) - n(20)), &
                   "n axis must be compressed towards half filling")
        call check(n(2) - n(1) < n(21) - n(20), &
                   "n axis must be compressed at the low-density end")
        ! Both ends are compressed, but not equally: `grade_high` (1.45) is the
        ! larger exponent, so half filling must end up finer than the dilute
        ! end.  Without this the two exponents could be swapped unnoticed.
        call check(n(40) - n(39) < 0.5_dp * (n(2) - n(1)), &
                   "the n axis must be finer at half filling than at n_min, " // &
                   "i.e. n_grade_high > n_grade_low")
    end subroutine test_graded_grids

    !> Regression test for the two round-off cliffs of `compute_V_xc_numerical`.
    !!
    !! (a) `n = 1` is a **discontinuity** of `V_xc` (the Mott cusp).  The
    !! direct generator evaluator and the spline consumer share
    !! `HALF_FILLING_SNAP_TOL`: round-off inside that band stays on the lower
    !! limit, while a genuine point outside it uses particle-hole symmetry.
    !! Generated tables are restricted to `n <= 1`, so no table node can be
    !! evaluated by the two modules on opposite sides of this cusp.
    !!
    !! (b) On the other side, `n - m -> 0` shrinks the `n` stencil to
    !! `0.3 (n - m)`, whose quotient amplifies the solver tolerance without
    !! bound.  Inside `M_EDGE_TOL` the exact edge identity is used instead, so
    !! a point at `n - m = 1e-9` must come out with `V_up = 0` exactly, not
    !! with amplified noise.
    subroutine test_half_filling_and_edge_guards()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: compute_V_xc_numerical, xc_potentials_t
        use lsda_constants, only: dp, HALF_FILLING_SNAP_TOL

        type(xc_potentials_t) :: v_half, v_snap, v_over, v_under, v_edge, v_near_edge

        ! (a) half filling reached from above by round-off shares the lower
        ! one-sided limit.  Outside the common band, particle-hole symmetry
        ! must select the upper one-sided limit instead.
        v_half = compute_V_xc_numerical(0.5_dp, 0.5_dp, 4.0_dp)
        v_snap = compute_V_xc_numerical(0.5_dp + 0.10_dp * HALF_FILLING_SNAP_TOL, &
                                        0.5_dp + 0.10_dp * HALF_FILLING_SNAP_TOL, 4.0_dp)
        v_over = compute_V_xc_numerical(0.5_dp + 1.0e-9_dp, 0.5_dp + 1.0e-9_dp, &
                                        4.0_dp)
        v_under = compute_V_xc_numerical(0.5_dp - 1.0e-9_dp, 0.5_dp - 1.0e-9_dp, &
                                         4.0_dp)

        call check(abs(v_snap%v_xc_up - v_half%v_xc_up) < 1.0e-8_dp, &
                   "a round-off excursion inside the shared snap band must keep V_xc_up")
        call check(abs(v_snap%v_xc_down - v_half%v_xc_down) < 1.0e-8_dp, &
                   "a round-off excursion inside the shared snap band must keep V_xc_down")
        call check(abs(v_over%v_xc_up + v_under%v_xc_up) < 1.0e-4_dp, &
                   "outside the shared snap band V_xc_up must obey particle-hole symmetry")
        call check(abs(v_over%v_xc_down + v_under%v_xc_down) < 1.0e-4_dp, &
                   "outside the shared snap band V_xc_down must obey particle-hole symmetry")

        ! (b) approach to the fully polarized edge
        v_edge = compute_V_xc_numerical(0.5_dp, 0.0_dp, 4.0_dp)
        v_near_edge = compute_V_xc_numerical(0.5_dp, 1.0e-9_dp, 4.0_dp)

        call check(abs(v_near_edge%v_xc_up) < 1.0e-14_dp, &
                   "a point 1e-9 from the m = n edge must use the edge identity, " // &
                   "not a stencil that amplifies the solver tolerance")
        call check(abs(v_near_edge%v_xc_down - v_edge%v_xc_down) < 1.0e-6_dp, &
                   "V_xc_down must be continuous into the fully polarized edge")
    end subroutine test_half_filling_and_edge_guards

    !> The generator must cover the full interaction range of the integral solver.
    !!
    !! The smallest accepted magnitude, `|U| = 0.5`, has no external-table
    !! oracle.  It is therefore checked against an intentionally over-resolved
    !! `Lambda` quadrature on a row containing the small nonzero magnetization
    !! node.  The comparison is relative and on the exchange splitting (see
    !! `check_table_convergence`), which is the quantity the Lambda floor
    !! exists for.  Both signs are generated because the attractive branch
    !! reaches the same repulsive solver through the Shiba map.
    !!
    !! `U = 0` remains refused: a zero XC table has no consumer, whereas point
    !! evaluation at `U = 0` remains legal (`e_xc == 0`), as covered by
    !! `compute_E_xc_U0`.
    subroutine test_table_generator_u_floor()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: generate_xc_table, grid_params_t, U_TABLE_MIN
        use lsda_errors, only: ERROR_INVALID_INPUT
        use table_io, only: xc_table_t
        use lsda_constants, only: dp

        type(grid_params_t) :: params, fine_params
        type(xc_table_t) :: table, fine_table
        integer :: status

        call check(abs(U_TABLE_MIN - 0.5_dp) < TOL, &
                   "the table floor must match the integral solver at |U| = 0.5")

        params%n_min = 0.5_dp
        params%n_max = 0.5_dp
        params%n_points = 1
        params%m_points = 3
        fine_params = params
        ! Must exceed n_lambda_floor(U/4) = 40 at U = 0.5, otherwise
        ! n_lambda_per_panel clamps both runs to the same order and the
        ! comparison is vacuous.
        fine_params%quad%n_lambda = 64

        call generate_xc_table(0.5_dp, params, table, status)
        call check(status == 0, "the repulsive U = 0.5 table must generate")
        call generate_xc_table(0.5_dp, fine_params, fine_table, status)
        call check(status == 0, "the over-resolved repulsive U = 0.5 table must generate")
        call check(all(ieee_is_finite(table%exc)) .and. all(ieee_is_finite(table%vxc_up)) &
                   .and. all(ieee_is_finite(table%vxc_down)), &
                   "the repulsive U = 0.5 table must contain only finite values")
        call check_table_convergence(table, fine_table, "U = 0.5")
        deallocate(table%n_grid, table%m_grid, table%exc, table%vxc_up, table%vxc_down)
        deallocate(fine_table%n_grid, fine_table%m_grid, fine_table%exc, fine_table%vxc_up, fine_table%vxc_down)

        call generate_xc_table(-0.5_dp, params, table, status)
        call check(status == 0, "the attractive U = -0.5 table must generate")
        call generate_xc_table(-0.5_dp, fine_params, fine_table, status)
        call check(status == 0, "the over-resolved attractive U = -0.5 table must generate")
        call check(all(ieee_is_finite(table%exc)) .and. all(ieee_is_finite(table%vxc_up)) &
                   .and. all(ieee_is_finite(table%vxc_down)), &
                   "the attractive U = -0.5 table must contain only finite values")
        call check_table_convergence(table, fine_table, "U = -0.5")
        deallocate(table%n_grid, table%m_grid, table%exc, table%vxc_up, table%vxc_down)
        deallocate(fine_table%n_grid, fine_table%m_grid, fine_table%exc, fine_table%vxc_up, fine_table%vxc_down)

        call generate_xc_table(0.0_dp, params, table, status)
        call check(status == ERROR_INVALID_INPUT, &
                   "U = 0 must be refused for table generation as well, not " // &
                   "silently treated as a special case")
    end subroutine test_table_generator_u_floor

    !> Compare a production table against an over-resolved one, *relatively*.
    !!
    !! The quantity that decides the Stoner instability is the exchange
    !! splitting `V_up - V_dn`, and at the smallest `m` node it is of order
    !! 1e-6.  The previous absolute tolerance of 1e-6 on the potentials was
    !! therefore about a thousand times the signal: forcing the Lambda order
    !! down to 20 gave a 12% error in the splitting and the test still passed,
    !! it only failed once the sign flipped.  Comparing the splitting relatively
    !! - and requiring the sign to survive - is what the test is for.
    !!
    !! @param[in] table  Table produced with the production quadrature
    !! @param[in] fine   Table produced with a deliberately higher Lambda order
    !! @param[in] label  Text identifying the case in the failure message
    subroutine check_table_convergence(table, fine, label)
        use fortuno_serial, only: check => serial_check
        use table_io, only: xc_table_t
        use lsda_constants, only: dp

        type(xc_table_t), intent(in) :: table, fine
        character(len=*), intent(in) :: label

        real(dp) :: split_p, split_f, rel
        integer :: i, j

        call check(maxval(abs(table%exc - fine%exc)) < 1.0e-6_dp, &
                   "the production " // label // " e_xc must self-converge")

        do i = 1, table%n_points_n
            do j = 1, table%n_points_m
                split_p = table%vxc_up(j, i) - table%vxc_down(j, i)
                split_f = fine%vxc_up(j, i) - fine%vxc_down(j, i)
                ! The m = 0 node has no splitting at all; nothing to compare.
                if (abs(split_f) < 1.0e-14_dp) cycle
                rel = abs(split_p - split_f) / abs(split_f)
                call check(split_p * split_f > 0.0_dp, &
                           "the production " // label // " exchange splitting must " // &
                           "keep the sign of the over-resolved one")
                call check(rel < 0.02_dp, &
                           "the production " // label // " exchange splitting must " // &
                           "converge relatively against the over-resolved table")
            end do
        end do
    end subroutine check_table_convergence

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
        U = 2.0_dp
        
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
