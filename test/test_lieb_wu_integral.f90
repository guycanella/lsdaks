!> Unit tests for the thermodynamic-limit Lieb-Wu solver.
!!
!! The anchors are, in order of strength:
!! 1. The closed form of the half-filled unpolarized ground state,
!!    `e = -4 Int_0^inf J0(w) J1(w) / (w (1 + exp(w U / 2))) dw`, evaluated by
!!    an independent composite quadrature inside this file.
!! 2. The exactly solvable limits: `U -> 0` (free gas, `e_xc = 0`),
!!    `U -> inf` (spinless fermions, `e = -(2/pi) sin(pi n)`), full
!!    polarization (`e_xc = 0`).
program test_lieb_wu_integral
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp, PI
    use lieb_wu_integral, only: lw_quad_t, gauss_legendre, &
                                lieb_wu_energy, lieb_wu_exc, lieb_wu_solve_qb, &
                                lieb_wu_solve_unpolarized
    use bethe_tables, only: compute_V_xc_numerical, xc_potentials_t
    use lsda_errors, only: ERROR_SUCCESS, ERROR_CONVERGENCE_FAILED
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
            test("lambda_mesh_capacity_failure", test_lambda_mesh_capacity_failure) &
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
    !!
    !! `fine%n_lambda` has to exceed `n_lambda_floor`, which is 40 at `U = 0.5`:
    !! `n_lambda_per_panel` takes `max(nominal, floor)`, so a "fine" order below
    !! the floor gives both runs exactly the same `Lambda` panels and the
    !! comparison silently degenerates into a test of `n_k` and `n_omega` alone
    !! - at the one `U` where the `Lambda` quadrature is the delicate part.
    subroutine test_quadrature_convergence()
        use fortuno_serial, only: check => serial_check

        type(lw_quad_t) :: quad, fine
        real(dp) :: exc_default, exc_fine
        integer :: ierr

        fine%n_k = 96
        fine%n_lambda = 64
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
    !! motivated it: the exchange splitting `V_up - V_down` near `m = 0`.
    !!
    !! This sweeps `U` instead of testing a single value, because the order the
    !! spin quadrature needs is not monotone in `U` (see `n_lambda_floor`): the
    !! leftover dyadic panel at the end of the `Lambda` mesh makes it oscillate,
    !! so a ladder calibrated on a handful of `U` can be badly wrong *between*
    !! its calibration points.  The grid deliberately includes the interior of
    !! every branch of the floor - `1.02` in particular, where the previous
    !! ladder returned the splitting with the **wrong sign** - and the hard band
    !! around `U = 0.65` to `0.9`.
    !!
    !! The comparison is relative and includes the sign, because the splitting
    !! at `m/n = 1e-5` is of order 1e-6 and any absolute tolerance large enough
    !! to be comfortable is larger than the signal itself.
    subroutine test_low_m_quadrature_floor_convergence()
        use fortuno_serial, only: check => serial_check

        type(lw_quad_t) :: fine
        type(xc_potentials_t) :: vp, vf
        real(dp), parameter :: U_LIST(9) = [0.5_dp, 0.65_dp, 0.9_dp, 1.0_dp, 1.02_dp, &
                                            1.1_dp, 1.3_dp, 1.5_dp, 2.0_dp]
        real(dp), parameter :: N_LIST(3) = [0.3_dp, 0.8_dp, 0.95_dp]
        real(dp), parameter :: M_FRAC(2) = [1.0e-5_dp, 1.0e-4_dp]
        real(dp) :: m, n, split_p, split_f, rel_err
        integer :: iu, in_, im

        ! Every quadrature dimension is refined, not just `Lambda`: with `n_k`
        ! and `n_omega` left at their production defaults the reference shares
        ! the charge and Fourier quadrature with the value under test, so the
        ! comparison is blind to a deficiency there and can only see `Lambda`.
        fine%n_lambda = 64
        fine%n_k = 96
        fine%n_omega = 64
        do iu = 1, size(U_LIST)
            do in_ = 1, size(N_LIST)
                do im = 1, size(M_FRAC)
                    n = N_LIST(in_)
                    m = n * M_FRAC(im)
                    vp = compute_V_xc_numerical(0.5_dp * (n + m), 0.5_dp * (n - m), &
                                                U_LIST(iu), delta_n=0.3_dp * m)
                    vf = compute_V_xc_numerical(0.5_dp * (n + m), 0.5_dp * (n - m), &
                                                U_LIST(iu), quad=fine, delta_n=0.3_dp * m)
                    split_p = vp%v_xc_up - vp%v_xc_down
                    split_f = vf%v_xc_up - vf%v_xc_down
                    rel_err = abs(split_p - split_f) / max(abs(split_f), 1.0e-12_dp)
                    call check(split_p * split_f > 0.0_dp, &
                               "production low-m splitting must preserve the fine-grid sign")
                    call check(rel_err < 0.02_dp, &
                               "production low-m splitting must converge against n_lambda=64")
                end do
            end do
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

    !> An extreme finite spin cutoff must return an error, not terminate the process.
    subroutine test_lambda_mesh_capacity_failure()
        use fortuno_serial, only: check => serial_check

        type(lw_quad_t) :: quad
        real(dp) :: n_tot, n_down, e_ba
        integer :: ierr

        call lieb_wu_solve_qb(1.0_dp, 1.0e6_dp, 0.5_dp, quad, n_tot, n_down, e_ba, ierr)
        call check(ierr == ERROR_CONVERGENCE_FAILED, &
                   "an over-capacity Lambda mesh must report convergence failure")
    end subroutine test_lambda_mesh_capacity_failure

end program test_lieb_wu_integral
