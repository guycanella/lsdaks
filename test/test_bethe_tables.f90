!> Unit tests for bethe_tables module
program test_bethe_tables
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
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
            test("compute_V_xc_polarized_corner_u4", test_compute_V_xc_polarized_corner_u4), &
            test("attractive_shiba_internal_consistency_u4", test_attractive_shiba_u4), &
            test("vxc_edge_identity_only_on_edge_u4", test_vxc_edge_identity_only_on_edge), &
            test("graded_grids", test_graded_grids), &
            test("low_m_exchange_splitting", test_low_m_exchange_splitting), &
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

    !> Regression test for the fully polarized corner of the U=4 table.
    !!
    !! At `n = 1`, `m = 1` the only stencils that stay inside the physical
    !! triangle are the two edges: `e_xc == 0` along `m = n` forces
    !! `V_up = 0`, and the `m` derivative along `n = 1` gives
    !! `V_dn = -2 de/dm = -4 (sqrt(2) - 1)` for U = 4.
    subroutine test_compute_V_xc_polarized_corner_u4()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: compute_E_xc, compute_V_xc_numerical, xc_potentials_t
        use lsda_constants, only: dp

        real(dp), parameter :: CPP_VXC_DN = -1.656854249466_dp
        type(xc_potentials_t) :: v_xc
        real(dp) :: e_xc

        e_xc = compute_E_xc(1.0_dp, 0.0_dp, 4.0_dp)
        v_xc = compute_V_xc_numerical(1.0_dp, 0.0_dp, 4.0_dp)

        call check(abs(e_xc) < 1.0e-14_dp, "e_xc must vanish at the polarized corner")
        call check(abs(v_xc%v_xc_up) < 1.0e-14_dp, &
                   "V_xc_up must vanish at the polarized corner")
        call check(abs(v_xc%v_xc_down - CPP_VXC_DN) < 1.0e-6_dp, &
                   "V_xc_down at the polarized corner must match the C++ reference")
    end subroutine test_compute_V_xc_polarized_corner_u4

    !> **Internal consistency only** of the attractive Shiba branch.
    !!
    !! `e_xc(n_up, n_dn; -U) = e_xc(1 - n_up, n_dn; U)`, with `V_xc_up`
    !! changing sign and `V_xc_down` keeping it.  Both sides of this test go
    !! through the same `shiba_map`, so it pins the convention the module uses
    !! internally and nothing more: flipping the channel that changes sign in
    !! **both** the generator and `xc_lsda` would leave this test passing.  The
    !! external anchor for the attractive convention is
    !! `test_lieb_wu_integral / reference_table_u_minus4`, which compares
    !! against the converted C++ table; do not weaken that one on the strength
    !! of this one.
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

    !> Regression test: the `m -> 0` exchange splitting at low `U`.
    !!
    !! @note VALIDATION-ONLY: this test and its helper `low_m_splitting` are
    !!       written against the C++ reference tables
    !!       (`data/tables/fortran_native/xc_table_u{1,2}.00.dat`) and must be
    !!       deleted together with `original/` and `data/tables/`.  What they
    !!       pin (the sign and magnitude of the `m -> 0` splitting, blocker B1)
    !!       then has no oracle left, so the guard that has to survive is the
    !!       `n_lambda_floor` of `lieb_wu_integral` itself, whose docstring
    !!       records the measured numbers.
    !!
    !! `V_up - V_dn = 2 de_xc/dm` vanishes linearly as `m -> 0`, and its sign
    !! is what decides the Stoner instability of the SCF, so the smallest `m`
    !! nodes of a table have to get it right.  Two defects used to hide there,
    !! and neither is visible to an absolute tolerance - the quantity itself is
    !! only `~1e-6` at these nodes - so this test is written against the
    !! *relative* deviation from the reference table:
    !!
    !! 1. the default `Lambda` quadrature order resolved the `O(K m^2)` signal
    !!    only for `U >= 2`, so at `U = 1` the difference was pure quadrature
    !!    noise: the splitting came out with the **wrong sign**, four times too
    !!    large.  Pinned by the `U = 1` sign check.
    !! 2. the `m` step was fixed at `delta_n = 1e-4`, i.e. far larger than `m`
    !!    itself, so the stencil `[|m - h|, m + h]` sampled `[h, h]` and
    !!    returned the curvature averaged over `[0, 2h]` instead of its value
    !!    at `m`.  That bias is what the 2.5% bound at `U = 2` pins: the fixed
    !!    step gives 4.5% there against 1.2% for the shrinking one.
    !!
    !! `U = 1` is the weakest interaction with a C++ reference table; `j = 2, 3`
    !! are the two smallest non-zero nodes of the reference `m` ladder (`m / n`
    !! around 1e-5 and 1e-4).  Rows are sampled every eighth to keep the test
    !! near a second.
    subroutine test_low_m_exchange_splitting()
        use fortuno_serial, only: check => serial_check
        use lsda_constants, only: dp

        real(dp) :: worst
        logical :: signs_agree

        call low_m_splitting(1.0_dp, 'data/tables/fortran_native/xc_table_u1.00.dat', &
                             worst, signs_agree)
        call check(signs_agree, &
                   "the m -> 0 exchange splitting must not come out with the " // &
                   "wrong sign at U = 1")
        call check(worst < 0.1_dp, &
                   "the m -> 0 exchange splitting must match the U=1 reference " // &
                   "to 10% relative")

        call low_m_splitting(2.0_dp, 'data/tables/fortran_native/xc_table_u2.00.dat', &
                             worst, signs_agree)
        call check(signs_agree, "the U=2 splitting must keep its sign as well")
        call check(worst < 0.025_dp, &
                   "the m -> 0 exchange splitting must match the U=2 reference " // &
                   "to 2.5% relative, which needs a stencil step that shrinks with m")
    end subroutine test_low_m_exchange_splitting

    !> Worst relative deviation of `V_up - V_dn` from a reference table over
    !! the two smallest non-zero `m` nodes of every eighth density row.
    !!
    !! @note VALIDATION-ONLY: helper of `low_m_exchange_splitting`; it reads a
    !!       C++ reference table and goes away with `original/` and
    !!       `data/tables/`.
    !!
    !! @param[in]  U            Hubbard interaction
    !! @param[in]  path         Reference table to compare against
    !! @param[out] worst        Largest `|split / split_ref - 1|`
    !! @param[out] signs_agree  `.false.` if any node disagrees in sign
    subroutine low_m_splitting(U, path, worst, signs_agree)
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: compute_V_xc_numerical, xc_potentials_t
        use table_io, only: xc_table_t, read_fortran_table, deallocate_table
        use lsda_constants, only: dp

        real(dp), intent(in) :: U
        character(len=*), intent(in) :: path
        real(dp), intent(out) :: worst
        logical, intent(out) :: signs_agree

        type(xc_table_t) :: ref
        type(xc_potentials_t) :: v
        real(dp) :: n, m, split, split_ref
        integer :: ierr, i, j, n_eval

        worst = 0.0_dp
        signs_agree = .true.
        n_eval = 0

        call read_fortran_table(path, ref, ierr)
        call check(ierr == 0, "the low-m reference table must be readable")
        if (ierr /= 0) return

        do i = 1, ref%n_points_n, 8
            n = ref%n_grid(i)
            if (n <= 0.1_dp) cycle
            do j = 2, 3
                m = ref%m_grid(j, i)
                split_ref = ref%vxc_up(j, i) - ref%vxc_down(j, i)
                if (abs(split_ref) < 1.0e-300_dp) cycle
                ! Below the unpolarized cut-off of the generator the splitting
                ! is returned as exactly zero by spin symmetry (`M_SYMMETRY_TOL
                ! = 1e-6`), so those nodes have no relative deviation to test.
                if (m <= 1.0e-6_dp) cycle

                v = compute_V_xc_numerical(0.5_dp * (n + m), 0.5_dp * (n - m), U)
                split = v%v_xc_up - v%v_xc_down
                n_eval = n_eval + 1

                if (split * split_ref <= 0.0_dp) signs_agree = .false.
                worst = max(worst, abs(split / split_ref - 1.0_dp))
            end do
        end do

        call check(n_eval >= 10, "the low-m comparison must cover several rows")
        call deallocate_table(ref)
    end subroutine low_m_splitting

    !> Regression test for the two round-off cliffs of `compute_V_xc_numerical`.
    !!
    !! (a) `n = 1` is a **discontinuity** of `V_xc` (the Mott cusp), and table
    !! grids reach half filling through arithmetic that lands on
    !! `n = 1 + 2e-9`: the reference `U = 2` table has `n_grid(50) =
    !! 1.000000002`.  A bare `n > 1.0_dp` test sends that point through the
    !! particle-hole branch, which returns the *opposite* one-sided limit and
    !! swaps `V_up` with `V_dn` - the same defect listed as Bug #1 of
    !! `xc_lsda` in CLAUDE.md.
    !!
    !! (b) On the other side, `n - m -> 0` shrinks the `n` stencil to
    !! `0.3 (n - m)`, whose quotient amplifies the solver tolerance without
    !! bound.  Inside `M_EDGE_TOL` the exact edge identity is used instead, so
    !! a point at `n - m = 1e-9` must come out with `V_up = 0` exactly, not
    !! with amplified noise.
    subroutine test_half_filling_and_edge_guards()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: compute_V_xc_numerical, xc_potentials_t
        use lsda_constants, only: dp

        real(dp), parameter :: CPP_VXC = -0.643363524_dp
        type(xc_potentials_t) :: v_half, v_over, v_edge, v_near_edge

        ! (a) half filling reached from above by round-off
        v_half = compute_V_xc_numerical(0.5_dp, 0.5_dp, 4.0_dp)
        v_over = compute_V_xc_numerical(0.5_dp + 1.0e-9_dp, 0.5_dp + 1.0e-9_dp, &
                                        4.0_dp)

        ! The residual difference is 3.9e-5, not zero: the cusp folds the first
        ! stencil point from n = 1 + 2e-9 to n = 1 - 2e-9, which offsets the
        ! backward stencil by 4e-9 / h.  The defect being guarded against
        ! returns +0.643 instead, i.e. it is off by 1.3.
        call check(abs(v_over%v_xc_up - CPP_VXC) < 1.0e-4_dp, &
                   "n = 1 + 2e-9 must stay on the n < 1 side of the Mott cusp")
        call check(abs(v_over%v_xc_up - v_half%v_xc_up) < 1.0e-4_dp, &
                   "a round-off excursion above n = 1 must not flip V_xc_up")
        call check(abs(v_over%v_xc_down - v_half%v_xc_down) < 1.0e-4_dp, &
                   "a round-off excursion above n = 1 must not swap the spins")

        ! (b) approach to the fully polarized edge
        v_edge = compute_V_xc_numerical(0.5_dp, 0.0_dp, 4.0_dp)
        v_near_edge = compute_V_xc_numerical(0.5_dp, 1.0e-9_dp, 4.0_dp)

        call check(abs(v_near_edge%v_xc_up) < 1.0e-14_dp, &
                   "a point 1e-9 from the m = n edge must use the edge identity, " // &
                   "not a stencil that amplifies the solver tolerance")
        call check(abs(v_near_edge%v_xc_down - v_edge%v_xc_down) < 1.0e-6_dp, &
                   "V_xc_down must be continuous into the fully polarized edge")
    end subroutine test_half_filling_and_edge_guards

    !> The generator must refuse interactions it was never validated at.
    !!
    !! `lieb_wu_integral` answers down to `|U| = 0.5`, but below `|U| = 1`
    !! there is no reference table left and the `m -> 0` splitting is only
    !! accurate to tens of percent, so writing a table there would be a silent
    !! loss of accuracy.
    !!
    !! `U = 0` is refused as well, and that is checked here explicitly: the
    !! floor test is a plain `abs(U) < U_TABLE_MIN` with no carve-out, so a
    !! `U = 0` table would otherwise come back as a bare
    !! `ERROR_INVALID_INPUT` while the docs and the CLI message promised it
    !! worked.  Point evaluation at `U = 0` stays legal (`e_xc == 0`), which
    !! `compute_E_xc_U0` covers; only generation refuses.
    subroutine test_table_generator_u_floor()
        use fortuno_serial, only: check => serial_check
        use bethe_tables, only: generate_xc_table, grid_params_t, U_TABLE_MIN
        use lsda_errors, only: ERROR_INVALID_INPUT
        use table_io, only: xc_table_t
        use lsda_constants, only: dp

        type(grid_params_t) :: params
        type(xc_table_t) :: table
        integer :: status

        call check(abs(U_TABLE_MIN - 1.0_dp) < TOL, &
                   "the validated table floor is |U| = 1")

        params%n_min = 0.5_dp
        params%n_max = 0.5_dp
        params%n_points = 1
        params%m_points = 1

        call generate_xc_table(0.5_dp, params, table, status)
        call check(status == ERROR_INVALID_INPUT, &
                   "a table below the validated U floor must be refused")

        call generate_xc_table(-0.5_dp, params, table, status)
        call check(status == ERROR_INVALID_INPUT, &
                   "the attractive branch must obey the same U floor")

        call generate_xc_table(0.0_dp, params, table, status)
        call check(status == ERROR_INVALID_INPUT, &
                   "U = 0 must be refused for table generation as well, not " // &
                   "silently treated as a special case")
    end subroutine test_table_generator_u_floor

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
