!> Unit and validation tests for the thermodynamic-limit Lieb-Wu solver.
!!
!! The anchors are, in order of strength:
!! 1. The closed form of the half-filled unpolarized ground state,
!!    `e = -4 Int_0^inf J0(w) J1(w) / (w (1 + exp(w U / 2))) dw`, evaluated by
!!    an independent composite quadrature inside this file.
!! 2. The exactly solvable limits: `U -> 0` (free gas, `e_xc = 0`),
!!    `U -> inf` (spinless fermions, `e = -(2/pi) sin(pi n)`), full
!!    polarization (`e_xc = 0`).
!! 3. Node-by-node agreement with the reference tables converted from the C++
!!    implementation, for `U = 4` and, through the Shiba transformation, for
!!    `U = -4`.
!!
!! @note VALIDATION-ONLY: anchor 3 (`reference_table_u4`,
!!       `reference_table_u_minus4` and their helper
!!       `compare_against_reference`) exists only to validate this generator
!!       against the C++ implementation, and is meant to be deleted together
!!       with `original/` and `data/tables/` once the generated tables have
!!       been accepted - it reads a reference file that will not exist.
!!       Anchors 1 and 2 are self-contained and stay.  Grep for
!!       `VALIDATION-ONLY` to find every such block.
!!
!! @note The reference tables print the `m = 0` row with about 9 significant
!!       digits, so the achievable agreement there is around 1e-9 absolute,
!!       not machine precision.  The closed-form anchors are therefore
!!       checked against the analytic integral, not against the file.
program test_lieb_wu_integral
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp, PI
    use lieb_wu_integral, only: lw_quad_t, lw_seed_t, gauss_legendre, &
                                lieb_wu_energy, lieb_wu_exc, lieb_wu_solve_qb, &
                                lieb_wu_solve_unpolarized
    use bethe_tables, only: compute_V_xc_numerical, xc_potentials_t
    use lsda_errors, only: ERROR_SUCCESS
    implicit none

    call execute_serial_cmd_app(get_tests())

contains

    function get_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests

        tests = test_list([ &
            test("half_filling_closed_form", test_half_filling_closed_form), &
            test("weak_coupling_limit", test_weak_coupling_limit), &
            test("density_input_roundoff", test_density_input_roundoff), &
            test("strong_coupling_spinless", test_strong_coupling_spinless), &
            test("full_polarization_zero_exc", test_full_polarization_zero_exc), &
            test("particle_hole_and_spin_symmetry", test_particle_hole_and_spin_symmetry), &
            test("quadrature_convergence", test_quadrature_convergence), &
            test("low_m_quadrature_floor_convergence", test_low_m_quadrature_floor_convergence), &
            test("unpolarized_matches_large_b", test_unpolarized_matches_large_b), &
            test("reference_table_u4", test_reference_table_u4), &
            test("reference_table_u_minus4", test_reference_table_u_minus4) &
        ])
    end function get_tests

    !> `e = -4 Int_0^inf J0(w) J1(w) / (w (1 + exp(w U / 2))) dw`
    !!
    !! Independent composite Gauss-Legendre quadrature: the integrand decays
    !! like `exp(-w U / 2)` times an oscillation of period `2 pi`, so panels of
    !! length below one with 16 nodes each resolve it to machine precision.
    function half_filling_closed_form(U) result(e)
        real(dp), intent(in) :: U
        real(dp) :: e

        integer, parameter :: NP = 16
        real(dp) :: xx(NP), ww(NP), w_max, a, b, s
        integer :: p, n_panels, j

        w_max = 120.0_dp / U + 40.0_dp
        n_panels = int(4.0_dp * w_max) + 1
        s = 0.0_dp
        do p = 1, n_panels
            a = w_max * real(p - 1, dp) / real(n_panels, dp)
            b = w_max * real(p, dp) / real(n_panels, dp)
            call gauss_legendre(NP, a, b, xx, ww)
            do j = 1, NP
                s = s + ww(j) * bessel_j0(xx(j)) * bessel_j1(xx(j)) &
                        / (xx(j) * (1.0_dp + exp(0.5_dp * xx(j) * U)))
            end do
        end do
        e = -4.0_dp * s
    end function half_filling_closed_form

    !> Half-filled unpolarized energy against the closed form for U = 2, 4, 8.
    subroutine test_half_filling_closed_form()
        use fortuno_serial, only: check => serial_check

        real(dp), parameter :: U_LIST(3) = [2.0_dp, 4.0_dp, 8.0_dp]
        type(lw_quad_t) :: quad
        real(dp) :: e_num, e_exact
        integer :: ierr, i
        character(len=64) :: msg

        do i = 1, 3
            call lieb_wu_energy(1.0_dp, 0.0_dp, U_LIST(i), quad, e_num, ierr)
            call check(ierr == ERROR_SUCCESS, "half filling must solve")
            e_exact = half_filling_closed_form(U_LIST(i))
            write(msg, '(A,F4.1,A)') "n=1, m=0 energy must match the closed form at U=", &
                U_LIST(i), " to 1e-8"
            call check(abs(e_num - e_exact) < 1.0e-8_dp, trim(msg))
        end do
    end subroutine test_half_filling_closed_form

    !> `e_xc -> 0` as `U -> 0`, exactly 0 at U = 0, rejected below the floor.
    !!
    !! The decay is monotone and roughly quadratic in `U` down to the lowest
    !! interaction the single-panel `k` quadrature can resolve.  Below that
    !! floor the solver must refuse rather than return a wrong number - a
    !! regression there is what silently corrupted tables would look like.
    subroutine test_weak_coupling_limit()
        use fortuno_serial, only: check => serial_check

        real(dp), parameter :: U_LIST(3) = [2.0_dp, 1.0_dp, 0.5_dp]
        type(lw_quad_t) :: quad
        real(dp) :: exc(3), exc_zero
        integer :: ierr, i

        do i = 1, 3
            call lieb_wu_exc(0.3_dp, 0.2_dp, U_LIST(i), quad, exc(i), ierr)
            call check(ierr == ERROR_SUCCESS, "weak coupling must solve")
        end do

        call check(abs(exc(3)) < abs(exc(2)) .and. abs(exc(2)) < abs(exc(1)), &
                   "|e_xc| must shrink monotonically as U -> 0")
        call check(abs(exc(3)) < 3.0e-3_dp, "e_xc must be small at U = 0.5")

        call lieb_wu_exc(0.3_dp, 0.2_dp, 0.0_dp, quad, exc_zero, ierr)
        call check(ierr == ERROR_SUCCESS, "U = 0 must be accepted")
        call check(abs(exc_zero) < 1.0e-300_dp, "e_xc must be exactly 0 at U = 0")

        call lieb_wu_exc(0.3_dp, 0.2_dp, 0.05_dp, quad, exc_zero, ierr)
        call check(ierr /= ERROR_SUCCESS, &
                   "U below the quadrature floor must be rejected, not guessed")
    end subroutine test_weak_coupling_limit

    !> Small positive overshoots from a graded n=1 endpoint are accepted by
    !! the integral solver and then handled by the particle-hole clamp.
    subroutine test_density_input_roundoff()
        use fortuno_serial, only: check => serial_check

        type(lw_quad_t) :: quad
        real(dp) :: exc
        integer :: ierr

        call lieb_wu_exc(1.0_dp + 2.0e-9_dp, 0.0_dp, 4.0_dp, quad, exc, ierr)
        call check(ierr == ERROR_SUCCESS, &
                   "small positive density-grid overshoot must be accepted")
        call check(abs(exc) < 1.0e-10_dp, &
                   "the overshoot must map to the fully polarized zero-XC edge")
    end subroutine test_density_input_roundoff

    !> `U -> infinity` maps the unpolarized state onto spinless fermions.
    subroutine test_strong_coupling_spinless()
        use fortuno_serial, only: check => serial_check

        type(lw_quad_t) :: quad
        real(dp) :: e_num, e_spinless
        integer :: ierr

        call lieb_wu_energy(0.5_dp, 0.0_dp, 4.0e3_dp, quad, e_num, ierr)
        call check(ierr == ERROR_SUCCESS, "strong coupling must solve")
        e_spinless = -(2.0_dp / PI) * sin(PI * 0.5_dp)
        call check(abs(e_num - e_spinless) < 1.0e-3_dp, &
                   "large U must approach the spinless fermion energy")

        call lieb_wu_energy(0.3_dp, 0.0_dp, 4.0e3_dp, quad, e_num, ierr)
        e_spinless = -(2.0_dp / PI) * sin(PI * 0.3_dp)
        call check(abs(e_num - e_spinless) < 1.0e-3_dp, &
                   "large U must approach the spinless fermion energy at n=0.3")
    end subroutine test_strong_coupling_spinless

    !> The fully polarized edge, approached through the general branch.
    !!
    !! The exact-edge values (`e_xc = 0`, `e = -(2/pi) sin(pi n)`) come out of
    !! early returns in `lieb_wu_exc` and `lieb_wu_invert`, so asserting them on
    !! the edge alone tests nothing but those two `if` statements.  What has
    !! content is the **continuity** of the general branch with them: with a
    !! small but finite minority density the coupled `(Q, B)` inversion runs in
    !! full, with `B -> 0`, and its answer has to converge onto the closed form
    !! linearly in `n_down`.  A wrong `B -> 0` limit of the spin kernels shows
    !! up here and nowhere else.
    subroutine test_full_polarization_zero_exc()
        use fortuno_serial, only: check => serial_check

        real(dp), parameter :: N_TOT = 0.6_dp
        type(lw_quad_t) :: quad
        real(dp) :: exc, e_num, e_edge, d_coarse, d_fine, s_coarse, s_fine
        integer :: ierr

        e_edge = -(2.0_dp / PI) * sin(PI * N_TOT)

        ! Exact edge: the analytic shortcuts.
        call lieb_wu_exc(N_TOT, 0.0_dp, 4.0_dp, quad, exc, ierr)
        call check(ierr == ERROR_SUCCESS, "polarized state must solve")
        call check(abs(exc) < 1.0e-14_dp, "e_xc must vanish at n_down = 0")

        call lieb_wu_energy(N_TOT, N_TOT, 4.0_dp, quad, e_num, ierr)
        call check(abs(e_num - e_edge) < 1.0e-14_dp, &
                   "polarized energy must be the free spinless band")

        ! General branch at n_down = 1e-4 and 1e-6: B -> 0 without any shortcut.
        call lieb_wu_exc(N_TOT - 1.0e-4_dp, 1.0e-4_dp, 4.0_dp, quad, exc, ierr)
        call check(ierr == ERROR_SUCCESS, "near-polarized state must solve at 1e-4")
        call check(abs(exc) < 1.0e-3_dp, "e_xc must be small next to the edge")
        s_coarse = exc / 1.0e-4_dp
        call lieb_wu_energy(N_TOT, N_TOT - 2.0e-4_dp, 4.0_dp, quad, e_num, ierr)
        call check(ierr == ERROR_SUCCESS, "near-polarized energy must solve at 1e-4")
        d_coarse = e_num - e_edge

        call lieb_wu_exc(N_TOT - 1.0e-6_dp, 1.0e-6_dp, 4.0_dp, quad, exc, ierr)
        call check(ierr == ERROR_SUCCESS, "near-polarized state must solve at 1e-6")
        s_fine = exc / 1.0e-6_dp
        call lieb_wu_energy(N_TOT, N_TOT - 2.0e-6_dp, 4.0_dp, quad, e_num, ierr)
        call check(ierr == ERROR_SUCCESS, "near-polarized energy must solve at 1e-6")
        d_fine = e_num - e_edge

        call check(abs(d_coarse) < 2.0e-4_dp .and. abs(d_fine) < 2.0e-6_dp, &
                   "the general branch must converge onto the analytic edge")
        call check(abs(d_fine) < 0.02_dp * abs(d_coarse), &
                   "the approach to the edge must be linear in n_down")
        call check(abs(s_fine - s_coarse) < 1.0e-3_dp * abs(s_fine), &
                   "de_xc/dn_down must have a finite limit on the polarized edge")
    end subroutine test_full_polarization_zero_exc

    !> Particle-hole and spin-exchange symmetry of `e_xc`.
    !!
    !! `lieb_wu_exc` orders the two spin channels with `max`/`min` on its first
    !! two lines, so simply swapping `n_up` and `n_down` walks the identical
    !! code and can never fail.  The symmetry that does exercise a distinct
    !! path is the particle-hole one: `e_xc(n_up, n_dn) = e_xc(1 - n_up,
    !! 1 - n_dn)`.  The image point here has `n_up + n_dn = 1.3 > 1`, so it
    !! enters the reduction branch, gets mapped back into the triangle, and is
    !! solved at different rapidity cut-offs; a wrong particle-hole transform
    !! fails this test.  Both orderings of the image are checked so that the
    !! transform and the channel ordering are pinned together.
    subroutine test_particle_hole_and_spin_symmetry()
        use fortuno_serial, only: check => serial_check

        type(lw_quad_t) :: quad
        real(dp) :: exc_direct, exc_hole, exc_hole_swapped
        integer :: ierr

        call lieb_wu_exc(0.45_dp, 0.25_dp, 4.0_dp, quad, exc_direct, ierr)
        call check(ierr == ERROR_SUCCESS, "direct state must solve")
        call check(abs(exc_direct) > 1.0e-3_dp, &
                   "the reference point must have a non-trivial e_xc")

        call lieb_wu_exc(0.55_dp, 0.75_dp, 4.0_dp, quad, exc_hole, ierr)
        call check(ierr == ERROR_SUCCESS, "particle-hole image must solve")
        call check(abs(exc_hole - exc_direct) < 1.0e-12_dp, &
                   "e_xc must obey particle-hole symmetry above half filling")

        call lieb_wu_exc(0.75_dp, 0.55_dp, 4.0_dp, quad, exc_hole_swapped, ierr)
        call check(ierr == ERROR_SUCCESS, "swapped particle-hole image must solve")
        call check(abs(exc_hole_swapped - exc_direct) < 1.0e-12_dp, &
                   "e_xc must be invariant under particle-hole plus spin exchange")
    end subroutine test_particle_hole_and_spin_symmetry

    !> The default quadrature orders must already be converged.
    subroutine test_quadrature_convergence()
        use fortuno_serial, only: check => serial_check

        type(lw_quad_t) :: quad, fine
        real(dp) :: exc_default, exc_fine
        integer :: ierr

        fine%n_k = 96
        fine%n_lambda = 28
        fine%n_omega = 32

        call lieb_wu_exc(0.35_dp, 0.25_dp, 4.0_dp, quad, exc_default, ierr)
        call check(ierr == ERROR_SUCCESS, "default quadrature must solve")
        call lieb_wu_exc(0.35_dp, 0.25_dp, 4.0_dp, fine, exc_fine, ierr)
        call check(ierr == ERROR_SUCCESS, "fine quadrature must solve")
        call check(abs(exc_default - exc_fine) < 1.0e-12_dp, &
                   "default quadrature must be converged to 1e-12")

        call lieb_wu_exc(0.35_dp, 0.25_dp, 0.5_dp, quad, exc_default, ierr)
        call check(ierr == ERROR_SUCCESS, "default quadrature must solve at U = 0.5")
        call lieb_wu_exc(0.35_dp, 0.25_dp, 0.5_dp, fine, exc_fine, ierr)
        call check(ierr == ERROR_SUCCESS, "fine quadrature must solve at U = 0.5")
        call check(abs(exc_default - exc_fine) < 1.0e-10_dp, &
                   "default quadrature must be converged at the lowest allowed U")
    end subroutine test_quadrature_convergence

    !> The production Lambda-order floor must converge the observable that
    !! motivated it: the exchange splitting V_up - V_down near m=0. Comparing
    !! the default path with an intentionally over-resolved n_lambda=40 path
    !! keeps this regression independent of the C++ reference tables.
    subroutine test_low_m_quadrature_floor_convergence()
        use fortuno_serial, only: check => serial_check

        type(lw_quad_t) :: fine
        type(xc_potentials_t) :: vp, vf
        real(dp), parameter :: U = 1.0_dp, N = 0.5_dp
        real(dp), parameter :: M_FRAC(3) = [1.0e-5_dp, 1.0e-4_dp, 1.0e-3_dp]
        real(dp) :: m, split_p, split_f, rel_err
        integer :: i

        fine%n_lambda = 40
        do i = 1, size(M_FRAC)
            m = N * M_FRAC(i)
            vp = compute_V_xc_numerical(0.5_dp * (N + m), 0.5_dp * (N - m), U, &
                                         delta_n=0.3_dp * m)
            vf = compute_V_xc_numerical(0.5_dp * (N + m), 0.5_dp * (N - m), U, &
                                         quad=fine, delta_n=0.3_dp * m)
            split_p = vp%v_xc_up - vp%v_xc_down
            split_f = vf%v_xc_up - vf%v_xc_down
            rel_err = abs(split_p - split_f) / max(abs(split_f), 1.0e-12_dp)
            call check(split_p * split_f > 0.0_dp, &
                       "production low-m splitting must preserve the fine-grid sign")
            call check(rel_err < 0.02_dp, &
                       "production low-m splitting must converge against n_lambda=40")
        end do
    end subroutine test_low_m_quadrature_floor_convergence

    !> The `B = infinity` Shiba branch must agree with a large finite `B`.
    !!
    !! The magnetization vanishes as `exp(-pi B / (2 u))`, so at `B = 20` with
    !! `u = 1` the finite-B solve is unpolarized to far below the working
    !! precision and both branches must return the same energy.  This is the
    !! check that the Fourier kernel `R(x)` of the reduction is right.
    subroutine test_unpolarized_matches_large_b()
        use fortuno_serial, only: check => serial_check

        type(lw_quad_t) :: quad
        real(dp) :: n_inf, e_inf, n_fin, ndn_fin, e_fin
        integer :: ierr

        call lieb_wu_solve_unpolarized(1.9_dp, 4.0_dp, quad, n_inf, e_inf, ierr)
        call check(ierr == ERROR_SUCCESS, "unpolarized branch must solve")
        call lieb_wu_solve_qb(1.9_dp, 20.0_dp, 4.0_dp, quad, n_fin, ndn_fin, e_fin, ierr)
        call check(ierr == ERROR_SUCCESS, "finite-B branch must solve")

        call check(abs(n_inf - n_fin) < 1.0e-12_dp, &
                   "density of both branches must agree at large B")
        call check(abs(e_inf - e_fin) < 1.0e-12_dp, &
                   "energy of both branches must agree at large B")
        call check(abs(n_fin - 2.0_dp * ndn_fin) < 1.0e-12_dp, &
                   "large B must be unpolarized")
    end subroutine test_unpolarized_matches_large_b

    !> Compare `exc`, `Vxc_up` and `Vxc_down` on the nodes of the U = 4 reference.
    !!
    !! @note VALIDATION-ONLY: delete this test together with `original/` and
    !!       `data/tables/`; it reads
    !!       `data/tables/fortran_native/xc_table_u4.00.dat`.
    !!
    !! The reference grid is non-uniform on both axes, so the comparison is
    !! made at the nodes of the file itself rather than on a grid of our own.
    !! The `m = 0` and `m = n` rows are always included: they are the two
    !! boundary identities (`de/dm = 0` and `V_up = 0`) that a wrong derivative
    !! stencil would break first.
    subroutine test_reference_table_u4()
        use fortuno_serial, only: check => serial_check

        real(dp) :: worst_exc, worst_up, worst_dn, worst_split
        character(len=96) :: split_msg

        call compare_against_reference(4.0_dp, worst_exc, worst_up, worst_dn, worst_split)

        call check(worst_exc < 1.0e-6_dp, "U=4 exc must match the reference to 1e-6")
        call check(worst_up < 1.0e-4_dp, "U=4 Vxc_up must match the reference to 1e-4")
        call check(worst_dn < 1.0e-4_dp, "U=4 Vxc_down must match the reference to 1e-4")
        write(split_msg, '(A,F10.4)') "U=4 exchange splitting relative error = ", worst_split
        call check(worst_split < 0.10_dp, trim(split_msg))
    end subroutine test_reference_table_u4

    !> Same comparison for `U = -4` through the Shiba transformation.
    !!
    !! @note VALIDATION-ONLY: delete this test together with `original/` and
    !!       `data/tables/`; it reads the same reference file.  The internal
    !!       Shiba consistency check
    !!       `test_bethe_tables / attractive_shiba_internal_consistency_u4`
    !!       survives the deletion but is *not* a replacement: it pins the
    !!       convention the module uses against itself, not against C++.
    !!
    !! The reference tables are stored for `|U|`, so the attractive values are
    !! defined by `exc(n_up, n_dn; -U) = exc(1 - n_up, n_dn; U)`,
    !! `Vxc_up = -Vxc_up(1 - n_up, n_dn; U)` and
    !! `Vxc_down = +Vxc_down(1 - n_up, n_dn; U)`.  The test walks the same
    !! reference nodes and feeds the generator the physical attractive
    !! densities that map onto them, so a sign slip in the mapping fails here.
    subroutine test_reference_table_u_minus4()
        use fortuno_serial, only: check => serial_check

        real(dp) :: worst_exc, worst_up, worst_dn

        call compare_against_reference(-4.0_dp, worst_exc, worst_up, worst_dn)

        call check(worst_exc < 1.0e-6_dp, "U=-4 exc must match the Shiba reference to 1e-6")
        call check(worst_up < 1.0e-4_dp, "U=-4 Vxc_up must match the Shiba reference to 1e-4")
        call check(worst_dn < 1.0e-4_dp, "U=-4 Vxc_down must match the Shiba reference to 1e-4")
    end subroutine test_reference_table_u_minus4

    !> Worst node-by-node deviation from the reference table of `|U|`.
    !!
    !! @note VALIDATION-ONLY: helper of the two reference-table tests above;
    !!       it goes away with them, `original/` and `data/tables/`.
    !!
    !! The reference `m` axis is graded: it walks five decades of `m / n` with
    !! only a handful of nodes per decade at the bottom and most of its nodes
    !! above `m = 0.1 n`.  A plain "every 20th node" stride therefore skips
    !! almost the whole refined region next to `m = 0`, which is exactly where
    !! `de_xc/dm` carries its logarithmic corrections and where a wrong
    !! derivative stencil would hide.  The sampled indices below are picked to
    !! land on the decade boundaries of that ladder (`j = 2, 4, 7, 12, 22`,
    !! i.e. `m / n` around 1e-5, 1.8e-4, 1e-3, 2e-3, 1e-2), on the interior,
    !! and on the last nodes before the fully polarized edge.  The full
    !! 10100-node comparison gives the same worst-case figures but takes
    !! minutes.
    !!
    !! @param[in]  U          Hubbard interaction, with sign
    !! @param[out] worst_exc  Largest `|exc - reference|`
    !! @param[out] worst_up   Largest `|Vxc_up - reference|`
    !! @param[out] worst_dn   Largest `|Vxc_down - reference|`
    !! @param[out] worst_split Optional normalized error in `Vxc_up - Vxc_down`
    !!                         near `m=0`, with a 1e-5 scale floor.
    subroutine compare_against_reference(U, worst_exc, worst_up, worst_dn, worst_split)
        use fortuno_serial, only: check => serial_check
        use table_io, only: xc_table_t, read_fortran_table, deallocate_table
        use bethe_tables, only: compute_E_xc, compute_V_xc_numerical, xc_potentials_t
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        real(dp), intent(in) :: U
        real(dp), intent(out) :: worst_exc, worst_up, worst_dn
        real(dp), intent(out), optional :: worst_split

        integer, parameter :: M_SAMPLE(13) = &
            [1, 2, 4, 7, 12, 22, 42, 82, 122, 162, 192, 199, 202]
        type(xc_table_t) :: ref
        type(xc_potentials_t) :: v
        type(lw_seed_t) :: seed
        real(dp) :: n, m, nu, nd, exc
        integer :: ierr, i, j, jj, n_eval
        logical :: all_finite

        worst_exc = 0.0_dp
        worst_up = 0.0_dp
        worst_dn = 0.0_dp
        if (present(worst_split)) worst_split = 0.0_dp

        call read_fortran_table('data/tables/fortran_native/xc_table_u4.00.dat', ref, ierr)
        call check(ierr == ERROR_SUCCESS, "U=4 reference table must be readable")
        if (ierr /= ERROR_SUCCESS) return

        all_finite = .true.
        n_eval = 0

        do i = 1, ref%n_points_n
            seed%valid = .false.
            do jj = 1, size(M_SAMPLE)
                j = min(M_SAMPLE(jj), ref%n_points_m)

                n = ref%n_grid(i)
                m = ref%m_grid(j, i)
                if (m > n) cycle

                ! Densities of the repulsive reference node; for U < 0 the
                ! Shiba map n_up -> 1 - n_up is what the generator undoes.
                nu = 0.5_dp * (n + m)
                nd = 0.5_dp * (n - m)
                if (U < 0.0_dp) nu = 1.0_dp - nu

                exc = compute_E_xc(nu, nd, U, seed=seed)
                v = compute_V_xc_numerical(nu, nd, U, seed=seed)
                n_eval = n_eval + 1

                if (.not. (ieee_is_finite(exc) .and. ieee_is_finite(v%v_xc_up) &
                           .and. ieee_is_finite(v%v_xc_down))) then
                    all_finite = .false.
                    cycle
                end if

                worst_exc = max(worst_exc, abs(exc - ref%exc(j, i)))
                if (U < 0.0_dp) then
                    worst_up = max(worst_up, abs(v%v_xc_up + ref%vxc_up(j, i)))
                else
                    worst_up = max(worst_up, abs(v%v_xc_up - ref%vxc_up(j, i)))
                end if
                worst_dn = max(worst_dn, abs(v%v_xc_down - ref%vxc_down(j, i)))
                if (present(worst_split)) then
                    if (m > 1.0e-6_dp .and. m / max(n, 1.0e-12_dp) <= 2.0e-2_dp) then
                        worst_split = max(worst_split, &
                            abs((v%v_xc_up - v%v_xc_down) - &
                                (ref%vxc_up(j, i) - ref%vxc_down(j, i))) / &
                            max(abs(ref%vxc_up(j, i) - ref%vxc_down(j, i)), 1.0e-5_dp))
                    end if
                end if
            end do
        end do

        call check(all_finite, "every reference node must evaluate to a finite value")
        call check(n_eval > 400, "the comparison must cover the whole reference grid")

        call deallocate_table(ref)
    end subroutine compare_against_reference

end program test_lieb_wu_integral
