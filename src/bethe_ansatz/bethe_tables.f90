!> Generation of exchange-correlation functional tables for BALDA.
!!
!! The XC energy per site comes from the thermodynamic-limit Lieb-Wu integral
!! equations (`lieb_wu_integral`), not from a finite-L discrete Bethe Ansatz
!! solve.  The tables consumed by the SCF are thermodynamic-limit functionals
!! `e_xc(n, m; U)`, so there is no finite-size error left in the generated
!! values.  The validated baseline for the default 75 x 202 grid is about 49 s
!! at `U = 4` on 14 OpenMP threads (306 s single-threaded).  Weak coupling is
!! the expensive end because the `Lambda` quadrature floor rises as `U` falls;
!! see `B_FLOOR_FRAC` in `lieb_wu_integral`.  The `b_max` cutoff is chosen from
!! numerical precision rather than a fixed wall-time promise, so timings depend
!! on the requested U.
!!
!! What this module adds on top of the integral-equation core:
!! 1. The attractive branch `U < 0`, via the Shiba partial particle-hole
!!    transformation (the same mapping `xc_lsda` applies when reading tables).
!! 2. `V_xc_sigma = d e_xc / d n_sigma` by finite differences, taken in the
!!    `(n, m)` basis so that every stencil stays inside the physical triangle
!!    `0 <= m <= n <= 1`.
!! 3. Grid assembly and file output in the reference table layout.
!!
!! @note The `exc` column has the Hartree term `U n_up n_dn` and the
!!       Kohn-Sham kinetic energy `T_s` already subtracted; the SCF adds
!!       `U n_other` back separately.  Do not double count.
!! @see lieb_wu_integral, table_io, xc_lsda
module bethe_tables
    use lsda_constants, only: dp, PI, TWOPI, U_SMALL, HALF_FILLING_SNAP_TOL
    use lieb_wu_integral, only: lw_quad_t, lw_seed_t, lieb_wu_exc
    use table_io, only: xc_table_t, write_fortran_table, count_nonfinite_entries, &
                        xc_table_filename
    use lsda_errors, only: ERROR_SUCCESS, ERROR_NOT_A_NUMBER, ERROR_INVALID_INPUT
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan, ieee_is_finite
    implicit none
    private

    !> Grid layout and numerical settings of a generated table.
    !!
    !! The magnetization axis runs from `0` to `n` inclusive, matching the
    !! reference tables and the expectation of `spline2d` (the `m < 0` and
    !! `n > 1` regions are reconstructed there by symmetry).  Both axes are
    !! graded, see `magnetization_grid` and `density_grid`; the defaults are
    !! the ones measured to bring the post-spline deviation from the reference
    !! `U = 4` table down to 2.7e-7 in `e_xc` and 7.7e-5 in `V_xc`.
    type, public :: grid_params_t
        real(dp) :: n_min = 0.02_dp      !< Minimum density
        real(dp) :: n_max = 1.0_dp       !< Maximum density
        integer :: n_points = 75         !< Number of density points
        integer :: m_points = 202        !< Number of magnetization points per density
        real(dp) :: m_frac_min = 1.0e-5_dp !< Smallest non-zero `m/n` of the graded axis
        real(dp) :: m_grade = 0.40_dp    !< Grading exponent of the `m` axis (1 = uniform)
        real(dp) :: n_grade_low = 1.1_dp  !< Grading exponent of the `n` axis at `n_min`
        real(dp) :: n_grade_high = 1.45_dp !< Grading exponent of the `n` axis at `n_max`
        real(dp) :: delta_n = 1.0e-4_dp  !< Finite-difference step of V_xc
        type(lw_quad_t) :: quad          !< Quadrature orders and inversion tolerances
    end type grid_params_t

    type, public :: xc_potentials_t
        real(dp) :: v_xc_up, v_xc_down
    end type xc_potentials_t

    !> Magnetizations below this are treated as exactly unpolarized, where
    !! `d e_xc / d m = 0` holds by spin symmetry.
    real(dp), parameter :: M_SYMMETRY_TOL = 1.0e-6_dp

    !> Relative distance from `m = n` within which the fully polarized edge
    !! identity `de/dn = -de/dm` is used instead of a stencil in `n`.
    !!
    !! The identity is exact only on the edge itself and carries an `O(n - m)`
    !! error away from it, so this must never be tied to the stencil step
    !! `delta_n` (doing so made interior points inherit `V_xc_up = 0`).  It is
    !! not pure round-off either: the `n` stencil has to shrink to
    !! `0.3 (n - m)` to stay inside the triangle, and its quotient amplifies
    !! the solver tolerance `quad%tol` by `1 / (0.6 (n - m))`.  The two error
    !! sources cross at `n - m ~ sqrt(tol) ~ 1e-6`, so below `1e-7` the
    !! identity is the more accurate of the two by an order of magnitude and
    !! the quotient is never evaluated there.  Every generated grid keeps its
    !! penultimate node at `n - m = 0.012 n`, four decades above this.
    real(dp), parameter :: M_EDGE_TOL = 1.0e-7_dp

    !> Relative width of the magnetization stencil next to `m = 0`.
    !!
    !! The step in `m` must shrink with `m`, or the stencil never samples the
    !! neighbourhood of `m` at all: with the fixed `h = 1e-4` it evaluated
    !! `[|m - h|, m + h] = [h, h]` and returned the curvature averaged over
    !! the stencil scale instead of the curvature at `m`, which is 2% wrong at
    !! `m / n = 1e-5` (measured against the reference `U = 4` table) and grows
    !! with the logarithmic corrections of the spin susceptibility as `U`
    !! drops.  `0.3 m` keeps the stencil inside the interval where `e_xc` is
    !! still quadratic in `m` while resolving `m` itself.
    real(dp), parameter :: M_STENCIL_FRAC = 0.3_dp

    !> Magnetization below which the stencil step stops shrinking with `m`.
    !!
    !! The signal the stencil has to resolve is
    !! `e(m + h_m) - e(m - h_m) = 4 K m h_m`, so a step proportional to `m`
    !! makes it fall as `m^2` and it disappears into the noise floor of
    !! `e_xc` - which is *not* the inversion tolerance (measured: tightening
    !! `quad%tol` from 1e-13 to 1e-15 changes nothing) but the discretization
    !! mismatch of the graded `Lambda` mesh between the two evaluations, whose
    !! panel count changes with `B ~ log(1/m)`.  Below this magnetization the
    !! step is therefore frozen at the width it has here (`0.3 M_STENCIL_MIN^2
    !! / m`, capped by the nominal `h`), which trades the noise for the much
    !! smaller bias of sampling `K` at the stencil scale.
    !!
    !! Measured worst relative error of `V_up - V_dn` over the smallest three
    !! `m` nodes of the reference tables, all `n` rows:
    !!
    !! | `U` | proportional step only | with this floor |
    !! |-----|------------------------|-----------------|
    !! |   1 |                  354%  |            8.2% |
    !! |   2 |                   84%  |            1.7% |
    !! |   4 |                 11.6%  |           0.92% |
    !! |   8 |                  1.4%  |           0.24% |
    real(dp), parameter :: M_STENCIL_MIN = 1.0e-5_dp

    !> Weakest interaction a *table* may be generated for.
    !!
    !! The point solver and generator share the lower limit
    !! `U_QUAD_MIN = 0.5`.  A table is more than a point evaluation: its
    !! smallest `m` nodes carry the `m -> 0` exchange splitting, whose
    !! logarithmic corrections sharpen as `U` falls.  The
    !! interaction-dependent `Lambda` quadrature floor in
    !! `lieb_wu_integral` resolves that splitting down to this limit; it is
    !! checked by self-convergence against an over-resolved quadrature rather
    !! than against an external table.
    !!
    !! `U = 0` is **not** carved out of this floor: `generate_xc_table` refuses
    !! it too.  `compute_E_xc` and `compute_V_xc_numerical` do accept `U = 0`
    !! and return exactly zero, but a whole *table* of zeros has no consumer -
    !! the SCF has no XC term at `U = 0` - so writing one would only look like
    !! a validated table without being one.
    real(dp), parameter, public :: U_TABLE_MIN = 0.5_dp

    public :: compute_E0
    public :: compute_E_xc
    public :: compute_V_xc_numerical
    public :: generate_xc_table
    public :: generate_table_grid
    public :: magnetization_grid
    public :: density_grid

contains

    !> Graded magnetization axis of one density row.
    !!
    !! A uniform `m` axis is not good enough for the splines that consume the
    !! table: the spin susceptibility of the gapless spin sector makes
    !! `de_xc/dm` carry logarithmic corrections as `m -> 0`, so a uniform axis
    !! leaves an O(1e-3) interpolation error in `V_xc_up - V_xc_dn` right next
    !! to `m = 0` - the region that dominates any Mott or antiferromagnetic
    !! calculation.  The reference tables converted from the C++ code grade
    !! their `m` axis over five decades for the same reason.
    !!
    !! The grading used here is a power law in the reduced variable `f = m / n`:
    !! the nodes are uniform in `s = f**grade`, so
    !!
    !! * `m(1) = 0` exactly (the unpolarized row, `de/dm = 0` by symmetry),
    !! * `m(2) = m_frac_min * n`,
    !! * `m(m_points) = n` exactly (the fully polarized edge),
    !! * the relative spacing `dm / m` scales as `m**(grade - 1)`.
    !!
    !! `grade = 1` is a uniform axis and `grade -> 0` a purely logarithmic one.
    !! The default `0.40` was chosen by measuring the post-spline deviation
    !! from the reference table: a purely logarithmic axis over-resolves
    !! `m < 1e-3 n`, where the functional is flat, and leaves the spacing near
    !! the `m = n` edge four times coarser than the reference, which showed up
    !! as an O(1e-3) error in `V_xc` at `m ~ 0.97 n`.
    !!
    !! @param[in]  n           Total density of the row, `> 0`
    !! @param[in]  m_frac_min  Smallest non-zero `m/n`, in `(0, 1)`
    !! @param[in]  grade       Grading exponent, in `(0, 1]`
    !! @param[out] m           Magnetization nodes; `size(m)` sets the count
    pure subroutine magnetization_grid(n, m_frac_min, grade, m)
        real(dp), intent(in) :: n, m_frac_min, grade
        real(dp), intent(out) :: m(:)

        integer :: j, m_points
        real(dp) :: s_min, s

        m_points = size(m)
        if (m_points < 1) return

        m(1) = 0.0_dp
        if (m_points == 1) return

        m(m_points) = n
        if (m_points == 2) return

        s_min = m_frac_min ** grade
        do j = 2, m_points - 1
            s = s_min + (1.0_dp - s_min) * real(j - 2, dp) / real(m_points - 2, dp)
            m(j) = n * s ** (1.0_dp / grade)
        end do
    end subroutine magnetization_grid

    !> Graded density axis of a table.
    !!
    !! A uniform `n` axis is not good enough either.  `e_xc(n)` has a cusp at
    !! half filling (the Mott gap opens there), which leaves an O(1e-4) error
    !! in `e_xc` and O(1e-3) in `V_xc` over the last two percent of the range;
    !! and `V_xc` turns over steeply as `n -> 0`, which leaves O(2e-4) in
    !! `V_xc` on the first few rows.  Both ends therefore need compression,
    !! which is what the reference tables do as well.  The nodes are placed by
    !! the monotone map
    !!
    !! `h(t) = t**a / (t**a + (1 - t)**b)`,  `t = (i - 1) / (n_points - 1)`,
    !!
    !! with `a = grade_low` controlling the `n_min` end and `b = grade_high`
    !! the `n_max` end; `a = b = 1` is the uniform axis and exponents above one
    !! push nodes towards the corresponding end.  The defaults `1.1` and `1.45`
    !! turn the spacing profile `0.020 / 0.023 / 0.020` of the uniform axis
    !! into `0.014 / 0.025 / 0.0035`; every figure here was measured as the
    !! post-spline deviation from the reference table, not guessed.
    !!
    !! @param[in]  n_min       First density node
    !! @param[in]  n_max       Last density node
    !! @param[in]  grade_low   Grading exponent at the `n_min` end, `>= 1`
    !! @param[in]  grade_high  Grading exponent at the `n_max` end, `>= 1`
    !! @param[out] n           Density nodes; `size(n)` sets the count
    pure subroutine density_grid(n_min, n_max, grade_low, grade_high, n)
        real(dp), intent(in) :: n_min, n_max, grade_low, grade_high
        real(dp), intent(out) :: n(:)

        integer :: i, n_points
        real(dp) :: t, h

        n_points = size(n)
        if (n_points < 1) return

        n(1) = n_min
        if (n_points == 1) return

        n(n_points) = n_max
        do i = 2, n_points - 1
            t = real(i - 1, dp) / real(n_points - 1, dp)
            h = t ** grade_low / (t ** grade_low + (1.0_dp - t) ** grade_high)
            n(i) = n_min + (n_max - n_min) * h
        end do
    end subroutine density_grid

    !> Non-interacting energy of a finite ring (U = 0, free Fermi gas).
    !!
    !! `E_0 = -2 sum_j cos(k_j)` for both spins with `k_j = (2 pi / L) I_j`,
    !! `I_j = j - (N+1)/2`.  This is the finite-L shell sum; the table
    !! generator uses the thermodynamic-limit form `-(2/pi) sin(pi n_sigma)`
    !! instead.  Kept for finite-size validation.
    !!
    !! @param[in] n_up  Spin-up density
    !! @param[in] n_dw  Spin-down density
    !! @param[in] L     System size
    !! @return          Non-interacting total energy E_0
    function compute_E0(n_up, n_dw, L) result(E0)
        real(dp), intent(in) :: n_up, n_dw
        integer, intent(in) :: L
        real(dp) :: E0, I_j, k_j
        integer :: Nup, Ndw, j

        Nup = NINT(n_up * real(L, dp))
        Ndw = NINT(n_dw * real(L, dp))

        E0 = 0.0_dp

        do j = 1, Nup
            I_j = real(j, dp) -0.5_dp * real(Nup + 1, dp)
            k_j = TWOPI * I_j / real(L, dp)
            E0 = E0 - 2.0_dp * cos(k_j)
        end do

        do j = 1, Ndw
            I_j = real(j, dp) -0.5_dp * real(Ndw + 1, dp)
            k_j = TWOPI * I_j / real(L, dp)
            E0 = E0 - 2.0_dp * cos(k_j)
        end do
    end function compute_E0

    !> XC energy per site, `e_xc = e_BA - T_s - U n_up n_dn`.
    !!
    !! For `U > 0` this delegates to the Lieb-Wu integral equations.  For
    !! `U < 0` the Shiba partial particle-hole transformation
    !! `e_xc(n_up, n_dn; U) = e_xc(1 - n_up, n_dn; |U|)` is applied first, so
    !! the attractive equations never have to be solved.
    !!
    !! @param[in]    n_up  Spin-up density in [0, 1]
    !! @param[in]    n_dn  Spin-down density in [0, 1]
    !! @param[in]    U     Hubbard interaction, with sign
    !! @param[in]    quad  Optional quadrature settings (defaults are used otherwise)
    !! @param[inout] seed  Optional warm start for the `(Q, B)` inversion
    !! @return             XC energy per site, NaN if the solve failed
    function compute_E_xc(n_up, n_dn, U, quad, seed) result(E_xc)
        real(dp), intent(in) :: n_up, n_dn, U
        type(lw_quad_t), intent(in), optional :: quad
        type(lw_seed_t), intent(inout), optional :: seed
        real(dp) :: E_xc

        type(lw_quad_t) :: q
        real(dp) :: nu, nd
        integer :: ierr

        E_xc = 0.0_dp
        if (abs(U) < U_SMALL) return

        if (present(quad)) q = quad

        call shiba_map(n_up, n_dn, U, nu, nd)

        if (present(seed)) then
            call lieb_wu_exc(nu, nd, abs(U), q, E_xc, ierr, seed)
        else
            call lieb_wu_exc(nu, nd, abs(U), q, E_xc, ierr)
        end if

        if (ierr /= ERROR_SUCCESS) E_xc = ieee_value(0.0_dp, ieee_quiet_nan)
    end function compute_E_xc

    !> XC potentials `V_xc_sigma = d e_xc / d n_sigma` by central differences.
    !!
    !! The derivatives are taken in the `(n, m)` basis and recombined as
    !! `V_up = de/dn + de/dm`, `V_dn = de/dn - de/dm`, because that is the
    !! basis in which every stencil can be kept inside the physical triangle
    !! `0 <= m <= n <= 1`:
    !!
    !! * `m = 0`: `de/dm = 0` by spin symmetry, no step in `m` needed.
    !! * `m = n` (fully polarized): `e_xc` vanishes identically along that
    !!   edge, hence `de/dn = -de/dm`, which reproduces `V_up = 0` and
    !!   `V_dn = -2 de/dm` of the reference tables.
    !! * `n = 1`: `e_xc` has a cusp at half filling (Mott gap), so the
    !!   derivative is the one-sided limit from `n < 1`, taken with a
    !!   second-order backward stencil.
    !!   Densities within `HALF_FILLING_SNAP_TOL` are snapped to that lower
    !!   limit. A genuine point above the band (for example `n = 1 + 2e-9`)
    !!   instead uses the particle-hole branch; for a polarized point, its
    !!   intermediate negative magnetization is spin-swapped when mapped back,
    !!   so the returned `V_xc_up`/`V_xc_down` are the upper-side channels.
    !!
    !! For `U < 0` the Shiba transformation is applied first and the spin-up
    !! potential changes sign, exactly as in `xc_lsda` and the C++ reference.
    !!
    !! @param[in]    n_up  Spin-up density in [0, 1]
    !! @param[in]    n_dn  Spin-down density in [0, 1]
    !! @param[in]    U     Hubbard interaction, with sign
    !! @param[in]    quad     Optional quadrature settings
    !! @param[inout] seed     Optional warm start for the `(Q, B)` inversion
    !! @param[in]    delta_n  Optional finite-difference step (default 1e-4)
    !! @return                XC potentials, NaN components if a solve failed
    function compute_V_xc_numerical(n_up, n_dn, U, quad, seed, delta_n) result(v_xc)
        real(dp), intent(in) :: n_up, n_dn, U
        type(lw_quad_t), intent(in), optional :: quad
        type(lw_seed_t), intent(inout), optional :: seed
        real(dp), intent(in), optional :: delta_n
        type(xc_potentials_t) :: v_xc

        type(lw_quad_t) :: q
        type(lw_seed_t) :: s
        real(dp) :: nu, nd, n, m, h, hm, hn, sn, sm, de_dn, de_dm, nan
        logical :: failed, on_edge

        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        v_xc%v_xc_up = 0.0_dp
        v_xc%v_xc_down = 0.0_dp
        if (abs(U) < U_SMALL) return

        if (present(quad)) q = quad
        if (present(seed)) s = seed

        call shiba_map(n_up, n_dn, U, nu, nd)

        ! Reduce to the canonical triangle, tracking how the derivatives of
        ! the canonical point map back: spin exchange flips de/dm, the full
        ! particle-hole transformation flips both de/dn and de/dm.
        n = nu + nd
        m = nu - nd
        sn = 1.0_dp
        sm = 1.0_dp
        ! Match the spline consumer at the Mott boundary: both the branch
        ! selection and the finite-difference stencil use the snapped density.
        ! Without this assignment, sharing the tolerance alone would still
        ! evaluate a near-boundary generator point at a different coordinate.
        if (abs(n - 1.0_dp) < HALF_FILLING_SNAP_TOL) n = 1.0_dp
        if (n > 1.0_dp + HALF_FILLING_SNAP_TOL) then
            n = 2.0_dp - n
            m = -m
            sn = -sn
            sm = -sm
        end if
        if (m < 0.0_dp) then
            m = -m
            sm = -sm
        end if

        if (n <= M_SYMMETRY_TOL) return

        h = 1.0e-4_dp
        if (present(delta_n)) h = delta_n
        h = min(h, 0.05_dp * n)
        failed = .false.
        on_edge = abs(n - m) <= M_EDGE_TOL * max(1.0_dp, n)

        ! --- d e_xc / d m at fixed n ---
        ! Shrink the step towards m = 0 as well, so that the stencil samples
        ! the neighbourhood of m instead of the fixed interval [0, 2h], and
        ! stop shrinking below M_STENCIL_MIN, where the difference would drown
        ! in the noise floor of e_xc.
        hm = min(h, M_STENCIL_FRAC &
                 * max(m, M_STENCIL_MIN**2 / max(m, tiny(1.0_dp))))
        if (m <= M_SYMMETRY_TOL) then
            de_dm = 0.0_dp
        else if (m + hm > n) then
            de_dm = (3.0_dp * exc_nm(n, m) - 4.0_dp * exc_nm(n, m - hm) &
                     + exc_nm(n, m - 2.0_dp * hm)) / (2.0_dp * hm)
        else
            de_dm = (exc_nm(n, m + hm) - exc_nm(n, abs(m - hm))) / (2.0_dp * hm)
        end if

        ! --- d e_xc / d n at fixed m ---
        if (on_edge) then
            ! e_xc == 0 all along m = n, so the directional derivative along
            ! that edge vanishes: de/dn + de/dm = 0.  This identity is exact
            ! ONLY on the edge itself, so it must be tied to the distance from
            ! the edge and never to the stencil step: at an interior point it
            ! would silently introduce an O(n - m) error.
            de_dn = -de_dm
        else
            ! Shrink the step so that the stencil cannot leave the physical
            ! triangle m <= n through the m = n edge.
            hn = min(h, 0.3_dp * (n - m))
            if (n + hn > 1.0_dp) then
                de_dn = (3.0_dp * exc_nm(n, m) - 4.0_dp * exc_nm(n - hn, m) &
                         + exc_nm(n - 2.0_dp * hn, m)) / (2.0_dp * hn)
            else
                de_dn = (exc_nm(n + hn, m) - exc_nm(n - hn, m)) / (2.0_dp * hn)
            end if
        end if

        if (failed) then
            v_xc%v_xc_up = nan
            v_xc%v_xc_down = nan
            return
        end if

        de_dn = sn * de_dn
        de_dm = sm * de_dm

        v_xc%v_xc_up = de_dn + de_dm
        v_xc%v_xc_down = de_dn - de_dm

        ! Shiba: the up potential of the attractive model is minus the up
        ! potential of the transformed repulsive one; the down potential keeps
        ! its sign.
        if (U < 0.0_dp) v_xc%v_xc_up = -v_xc%v_xc_up

        if (present(seed)) seed = s

    contains

        !> `e_xc` of the repulsive model at total density `nn` and magnetization `mm`.
        function exc_nm(nn, mm) result(e)
            real(dp), intent(in) :: nn, mm
            real(dp) :: e
            integer :: jerr

            call lieb_wu_exc(0.5_dp * (nn + mm), 0.5_dp * (nn - mm), abs(U), q, e, jerr, s)
            if (jerr /= ERROR_SUCCESS) then
                failed = .true.
                e = 0.0_dp
            end if
        end function exc_nm

    end function compute_V_xc_numerical

    !> Shiba partial particle-hole transformation for the attractive branch.
    !!
    !! `(n_up, n_dn; U < 0) -> (1 - n_up, n_dn; |U|)`.  For `U >= 0` the
    !! densities are passed through unchanged.
    !!
    !! @param[in]  n_up  Physical spin-up density
    !! @param[in]  n_dn  Physical spin-down density
    !! @param[in]  U     Hubbard interaction, with sign
    !! @param[out] nu    Transformed spin-up density
    !! @param[out] nd    Transformed spin-down density
    pure subroutine shiba_map(n_up, n_dn, U, nu, nd)
        real(dp), intent(in) :: n_up, n_dn, U
        real(dp), intent(out) :: nu, nd

        if (U < 0.0_dp) then
            nu = 1.0_dp - n_up
        else
            nu = n_up
        end if
        nd = n_dn
    end subroutine shiba_map

    !> Generate a complete XC table for one value of `U`.
    !!
    !! The density axis is the graded axis of `density_grid`; the magnetization axis
    !! of each row is the logarithmically graded axis of `magnetization_grid`,
    !! running from `0` to `n` inclusive.  The `m = n` end is the fully
    !! polarized state, where `e_xc` and `V_xc_up` vanish identically.  The
    !! converged `(Q, B)` pair of each point warm starts the next one, which is
    !! what keeps the whole grid within a minute of wall time (see the module
    !! header for the measured figures).
    !!
    !! Failed points are left as NaN; `write_fortran_table` refuses to persist
    !! a table that contains any.
    !!
    !! `|U| < U_TABLE_MIN` is refused, `U = 0` included: see `U_TABLE_MIN`.
    !!
    !! @param[in]  U      Hubbard interaction, with sign, `|U| >= U_TABLE_MIN`
    !! @param[in]  params Grid layout and numerical settings
    !! @param[out] table  Generated XC table
    !! @param[out] status `ERROR_SUCCESS`, or `ERROR_INVALID_INPUT` for any
    !!                    `|U| < U_TABLE_MIN` (`U = 0` too) or an inconsistent
    !!                    grid (the per-point failures show up as NaN)
    subroutine generate_xc_table(U, params, table, status)
        real(dp), intent(in) :: U
        type(grid_params_t), intent(in) :: params
        type(xc_table_t), intent(out) :: table
        integer, intent(out) :: status

        integer :: i, j
        real(dp) :: n, m, n_up, n_dw
        type(xc_potentials_t) :: v_xc
        type(lw_seed_t) :: seed

        status = ERROR_SUCCESS

        if (.not. ieee_is_finite(U) .or. abs(U) < U_TABLE_MIN) then
            status = ERROR_INVALID_INPUT
            return
        end if

        if (params%n_points < 1 .or. params%m_points < 1 .or. params%quad%n_k < 4 &
            .or. params%quad%n_lambda < 4 .or. params%quad%n_omega < 4 &
            .or. .not. ieee_is_finite(params%quad%tol) .or. params%quad%tol <= 0.0_dp &
            .or. .not. ieee_is_finite(params%delta_n) .or. params%delta_n <= 0.0_dp &
            .or. .not. ieee_is_finite(params%m_frac_min) .or. .not. ieee_is_finite(params%m_grade) &
            .or. .not. ieee_is_finite(params%n_grade_low) .or. .not. ieee_is_finite(params%n_grade_high) &
            .or. params%m_frac_min <= 0.0_dp .or. params%m_frac_min >= 1.0_dp &
            .or. params%m_grade <= 0.0_dp .or. params%m_grade > 1.0_dp &
            .or. params%n_grade_low < 1.0_dp .or. params%n_grade_high < 1.0_dp &
            .or. .not. ieee_is_finite(params%n_min) .or. .not. ieee_is_finite(params%n_max) &
            .or. params%n_min <= 0.0_dp .or. params%n_max > 1.0_dp &
            .or. params%n_min > params%n_max &
            .or. (params%n_points > 1 .and. params%n_min >= params%n_max)) then
            status = ERROR_INVALID_INPUT
            return
        end if

        table%U = U
        table%n_points_n = params%n_points
        table%n_points_m = params%m_points

        allocate(table%n_grid(params%n_points))
        allocate(table%m_grid(params%m_points, params%n_points))
        allocate(table%exc(params%m_points, params%n_points))
        allocate(table%vxc_up(params%m_points, params%n_points))
        allocate(table%vxc_down(params%m_points, params%n_points))

        table%exc = ieee_value(0.0_dp, ieee_quiet_nan)
        table%vxc_up = ieee_value(0.0_dp, ieee_quiet_nan)
        table%vxc_down = ieee_value(0.0_dp, ieee_quiet_nan)

        call density_grid(params%n_min, params%n_max, params%n_grade_low, &
                          params%n_grade_high, table%n_grid)

        !$OMP PARALLEL DO PRIVATE(i, j, n, m, n_up, n_dw, v_xc, seed) &
        !$OMP SHARED(table, params, U) SCHEDULE(dynamic)
        do i = 1, params%n_points
            n = table%n_grid(i)

            call magnetization_grid(n, params%m_frac_min, params%m_grade, &
                                    table%m_grid(:, i))

            seed = lw_seed_t()
            do j = 1, params%m_points
                m = table%m_grid(j, i)

                n_up = 0.5_dp * (n + m)
                n_dw = 0.5_dp * (n - m)

                table%exc(j, i) = compute_E_xc(n_up, n_dw, U, params%quad, seed)

                v_xc = compute_V_xc_numerical(n_up, n_dw, U, params%quad, seed, &
                                              params%delta_n)
                table%vxc_up(j, i) = v_xc%v_xc_up
                table%vxc_down(j, i) = v_xc%v_xc_down
            end do
        end do
        !$OMP END PARALLEL DO

    end subroutine generate_xc_table

    !> Generate and write XC tables for a list of `U` values.
    !!
    !! Each table is validated before it is written: a table with any
    !! non-finite entry aborts the sweep instead of being persisted, because
    !! the SCF splines would smear the invalid value over whole rows.
    !!
    !! @param[in]  U_values   Hubbard interaction values
    !! @param[in]  params     Grid layout and numerical settings
    !! @param[in]  output_dir Directory that receives the table files
    !! @param[out] status     `ERROR_SUCCESS` or the first failure code
    subroutine generate_table_grid(U_values, params, output_dir, status)
        real(dp), intent(in) :: U_values(:)
        type(grid_params_t), intent(in) :: params
        character(len=*), intent(in) :: output_dir
        integer, intent(out) :: status

        integer :: k, n_U, io_stat, n_bad_exc, n_bad_up, n_bad_dn
        real(dp) :: U_current
        type(xc_table_t) :: table
        character(len=256) :: filename

        status = ERROR_SUCCESS
        n_U = size(U_values)

        print '(A)', "========================================="
        print '(A)', "  XC Table Grid Generation"
        print '(A)', "========================================="
        print '(A,I0,A)', "Generating ", n_U, " tables"
        print '(A,F6.2,A,F6.2)', "U range: ", U_values(1), " to ", U_values(n_U)
        print '(A,I0,A,I0)', "Grid size: ", params%n_points, " x ", params%m_points
        print '(A)', ""

        do k = 1, n_U
            U_current = U_values(k)

            print '(A,I0,A,I0,A,F6.2)', "Processing U(", k, "/", n_U, ") = ", U_current

            call generate_xc_table(U_current, params, table, status)

            if (status /= ERROR_SUCCESS) then
                print '(A,F6.2)', "  ERROR: Failed to generate table for U = ", U_current
                return
            end if

            ! write_fortran_table refuses non-finite tables as well; the check
            ! is repeated here to name the offending quantity and stop the sweep
            ! before any later U is attempted.
            call count_nonfinite_entries(table, n_bad_exc, n_bad_up, n_bad_dn)
            if (n_bad_exc + n_bad_up + n_bad_dn > 0) then
                print '(A,F6.2,A)', "  ERROR: Table for U = ", U_current, &
                    " contains non-finite entries; not written."
                print '(A,I0,A,I0,A,I0)', "         non-finite: exc = ", n_bad_exc, &
                    ", Vxc_up = ", n_bad_up, ", Vxc_down = ", n_bad_dn
                status = ERROR_NOT_A_NUMBER
                return
            end if

            call xc_table_filename(trim(output_dir), U_current, filename)

            call write_fortran_table(filename, table, io_stat)

            if (io_stat /= ERROR_SUCCESS) then
                print '(A)', "  ERROR: Failed to write table file"
                status = io_stat
                return
            end if

            print '(A,A)', "  Saved: ", trim(filename)

            if (allocated(table%n_grid)) deallocate(table%n_grid)
            if (allocated(table%m_grid)) deallocate(table%m_grid)
            if (allocated(table%exc)) deallocate(table%exc)
            if (allocated(table%vxc_up)) deallocate(table%vxc_up)
            if (allocated(table%vxc_down)) deallocate(table%vxc_down)
        end do

        print '(A)', ""
        print '(A)', "========================================="
        print '(A,I0,A)', "Successfully generated ", n_U, " XC tables!"
        print '(A)', "========================================="

    end subroutine generate_table_grid
end module bethe_tables
