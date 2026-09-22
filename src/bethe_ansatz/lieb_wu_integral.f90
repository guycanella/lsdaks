!> Thermodynamic-limit Lieb-Wu integral equations for the 1D Hubbard model.
!!
!! This module replaces the finite-L discrete Bethe Ansatz solve as the
!! production path of the XC table generator.  Instead of solving for
!! `N_up + N_down` individual rapidities on a chain of `L` sites, it solves the
!! coupled integral equations for the rapidity *densities*, which are already
!! the `L -> infinity` limit.  The tabulated XC functional of BALDA is a
!! thermodynamic-limit quantity, so there is no `O(1/L)` or `O(1/L^2)` finite
!! size error left to fight.
!!
!! With `u = U/4` and `t = 1`, the ground state at density `n` and
!! magnetization `m = n_up - n_down` is described by a charge density `rho(k)`
!! on `[-Q, Q]` and a spin density `sigma(Lambda)` on `[-B, B]`:
!!
!! ```
!! rho(k)     = 1/(2 pi) + cos(k) * Int_{-B}^{B} a1(sin k - L') sigma(L') dL'
!! sigma(L)   = Int_{-Q}^{Q} a1(L - sin k) rho(k) dk
!!              - Int_{-B}^{B} a2(L - L') sigma(L') dL'
!! a_j(x)     = (1/2 pi) * 2 j u / ((j u)^2 + x^2)
!!
!! n          = Int_{-Q}^{Q} rho(k) dk
!! n_down     = Int_{-B}^{B} sigma(L) dL
!! e_BA       = -2 Int_{-Q}^{Q} cos(k) rho(k) dk
!! ```
!!
!! `e_BA` is the *total* ground state energy per site of
!! `H = -t sum c^+ c + U sum n_up n_down` (the interaction is encoded in the
!! charge rapidities), so the XC energy per site follows directly as
!!
!! ```
!! e_xc = e_BA - T_s - U n_up n_down,  T_s = -(2/pi) (sin(pi n_up) + sin(pi n_dn))
!! ```
!!
!! which is exactly the quantity stored in the `exc` column of the reference
!! tables (Hartree and Kohn-Sham kinetic energy already subtracted).
!!
!! Numerics
!! --------
!! * Both densities are even, so every integral is folded onto `[0, Q]` and
!!   `[0, B]` with symmetrized kernels.  This cuts the linear system size in
!!   half (and the cost by 8).
!! * `k` uses one Gauss-Legendre panel; `Lambda` uses dyadically graded
!!   Gauss-Legendre panels, because the kernel width is `O(u)` while `B` can
!!   reach `10^6` for `m -> 0`.  A single rule on `[0, B]` would not resolve
!!   the kernel there.
!! * `m = 0` means `B -> infinity`.  That limit is taken analytically with the
!!   Shiba/Lieb-Wu reduction: eliminating `sigma` by Fourier transform leaves a
!!   single equation in `k` with the kernel
!!   `R(x) = (1/2 pi) Int dw exp(i w x) / (1 + exp(2 u |w|))`.
!! * Full polarization (`n_down = 0`) is `B = 0`, `rho = 1/(2 pi)`, and
!!   `e_xc = 0` identically.
!! * The grid is specified in `(n, m)`, so `(Q, B)` is obtained by inversion.
!!   `n(Q, B)` is monotone in `Q` and `n_down(Q, B)` is monotone in `B`, so the
!!   inversion is a nested pair of bracketed root finds - it cannot run away.
!!   `B` is parameterized as `B = tan(w)`, `w` in `(0, pi/2)`, which is well
!!   scaled at both ends (`m -> 0` and `m -> n`).
!!
!! @note `U > 0` only.  Attractive interaction is handled by the caller with
!!       the Shiba partial particle-hole transformation (see `bethe_tables`
!!       and `xc_lsda`).
!! @note Validated range: `U >= 0.5`.  Below that the kernel width `u = U/4`
!!       shrinks faster than a fixed-order rule on `[0, Q]` can resolve, so
!!       `U < 0.5` is rejected with `ERROR_INVALID_INPUT` rather than
!!       answered with a silently wrong number.  Every tabulated `U` used by
!!       the SCF is at or above 1, and `|U| < U_SMALL` gives `e_xc = 0`
!!       exactly.
!! @see bethe_tables, xc_lsda
module lieb_wu_integral
    use lsda_constants, only: dp, PI, TWOPI, U_SMALL
    use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT, &
                           ERROR_CONVERGENCE_FAILED, ERROR_SINGULAR_MATRIX, &
                           ERROR_LAPACK_INVALID_ARG, ERROR_NOT_A_NUMBER
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    public :: lw_quad_t
    public :: lw_seed_t
    public :: gauss_legendre
    public :: lieb_wu_solve_qb
    public :: lieb_wu_solve_unpolarized
    public :: lieb_wu_invert
    public :: lieb_wu_energy
    public :: lieb_wu_exc

    !> Quadrature orders and root-finding tolerances.
    !!
    !! `n_k` is the node count on the *half* interval `[0, Q]`, i.e. the
    !! effective order on `[-Q, Q]` is `2 n_k`.  `n_lambda` is the node count
    !! per dyadic `Lambda` panel (typically 5 to 25 panels are used, so the
    !! effective order is again well above 64).  `n_omega` is the node count
    !! per unit panel of the Fourier integral of the `m = 0` kernel.
    type, public :: lw_quad_t
        integer  :: n_k = 40            !< Gauss-Legendre nodes on [0, Q]
        integer  :: n_lambda = 12       !< Gauss-Legendre nodes per Lambda panel
        integer  :: n_omega = 20        !< Gauss-Legendre nodes per omega panel
        real(dp) :: tol = 1.0e-13_dp    !< Density residual tolerance of the inversion
        integer  :: max_iter = 60       !< Maximum root-finder iterations per level
    end type lw_quad_t

    !> Warm start for the `(Q, w)` inversion, filled by `lieb_wu_invert`.
    !!
    !! Passing the converged pair of a neighbouring grid point cuts the number
    !! of linear solves per point by roughly a factor of three.
    type, public :: lw_seed_t
        logical  :: valid = .false.     !< `.true.` when `Q` and `w` are usable
        real(dp) :: Q = 0.0_dp          !< Charge rapidity cut-off
        real(dp) :: w = 0.0_dp          !< Spin cut-off angle, B = tan(w)
    end type lw_seed_t

    !> Below this magnetization the state is treated as unpolarized (B = inf).
    !! The energy error of the substitution is `O(m^2)`, i.e. below 1e-12.
    real(dp), parameter :: M_ZERO_TOL = 1.0e-6_dp
    !> Below this minority density the state is treated as fully polarized.
    real(dp), parameter :: N_MINORITY_TOL = 1.0e-14_dp
    !> Distance from `n = 1` within which `Q = pi` is used without inversion.
    real(dp), parameter :: N_FULL_TOL = 1.0e-13_dp
    !> Input-grid overshoot accepted before the caller's boundary clamp.
    !! This is wider than N_FULL_TOL because graded grids can produce a
    !! reproducible n = 1 + O(1e-9) endpoint through arithmetic roundoff.
    real(dp), parameter :: N_INPUT_TOL = 1.0e-8_dp
    !> Smallest admissible `B`; below it `sigma` is numerically zero.
    real(dp), parameter :: B_MIN = 1.0e-10_dp
    !> Maximum number of `Lambda` panels, including subdivisions of the last
    !! non-dyadic interval.
    integer, parameter :: MAX_LAMBDA_PANELS = 96
    !> Maximum width of a residual Lambda panel in units of `u`.
    !!
    !! The dyadic mesh can end at a point `v` for which `[v, B]` is several
    !! kernel widths wide. Splitting that residual prevents its geometry from
    !! changing the low-magnetization splitting as U varies.
    real(dp), parameter :: RESIDUAL_PANEL_MAX_WIDTH = 1.5_dp
    !> Exponential cut-off of the `m = 0` Fourier kernel: exp(-OMEGA_CUT).
    real(dp), parameter :: OMEGA_CUT = 48.0_dp
    !> Maximum number of `omega` panels of the `m = 0` Fourier kernel.
    integer, parameter :: MAX_OMEGA_PANELS = 96
    !> Smallest interaction the single-panel `k` quadrature can resolve.
    !!
    !! The kernels have width `u = U/4` in `sin k`, so as `U -> 0` they turn
    !! into delta functions that no fixed-order rule on `[0, Q]` can resolve.
    !! Measured against a 200-node reference, the default order is exact to
    !! 1e-12 at `U = 0.5`, to 5e-7 at `U = 0.2`, and diverges below `U = 0.05`.
    !! Rather than return silently wrong numbers, anything below this floor is
    !! rejected.  Every tabulated `U` of interest is at or above 1.
    real(dp), parameter, public :: U_QUAD_MIN = 0.5_dp

    !> Spin cut-off, in units of `u = U/4`, above which `n_lambda_floor` is
    !! applied; see `n_lambda_per_panel`.
    !!
    !! In the asymptotic regime `m / n ~ exp(-pi B / (2 u))`, so
    !! `B / u ~ (2/pi) ln(n/m)`: 2.94 is `m/n = 1e-2`, 4.4 is `1e-3`,
    !! and 7.33 is `1e-5`. At finite U the omitted prefactor involving
    !! `sin(Q)` shifts these values; the measured gate here is broader
    !! (`m/n` roughly 6e-3 to 8e-2 at U=1). Measured at `U = 1` on the
    !! default 75 x 202 grid (release, 14 OpenMP threads), against the worst
    !! relative error of the `m -> 0` splitting of the `U = 1` reference table
    !! (the blocker that `n_lambda_floor` exists for):
    !!
    !! | `B_FLOOR_FRAC` | table wall time | worst split error | sign flips |
    !! |----------------|-----------------|-------------------|------------|
    !! | 0 (always on)  |         106.8 s |             5.5%  |       0/75 |
    !! | 2.94           |          91.9 s |             5.5%  |       0/75 |
    !! | 4.4 (used)     |          81.6 s |             5.5%  |       0/75 |
    !! | 7.33           |          62.6 s |             5.5%  |       0/75 |
    !! | inf (never on) |          58.0 s |          9101   x |      30/75 |
    !!
    !! So the floor is indispensable, but only next to `m = 0`: switching it off
    !! above the measured gate costs nothing in the validation sweep and saves
    !! about 24% of the run. The chosen 4.4 is the smallest tested value that
    !! retains the required low-m splitting accuracy; 7.33 is not a safe default
    !! merely because it corresponds to a smaller asymptotic `m/n`.
    !!
    !! @note The wall times above were measured with the earlier, lower floor;
    !!       the current `n_lambda_floor` is more expensive below `U = 1.35`, so
    !!       treat them as the relative cost of the *gate*, not as absolute
    !!       timings. What the gate itself must satisfy - the low-`m` splitting
    !!       staying converged and sign-correct at `m/n` down to 1e-5 - is swept
    !!       over `U` and `n` by `test_low_m_quadrature_floor_convergence`.
    real(dp), parameter :: B_FLOOR_FRAC = 4.4_dp

    !> Residual of a scalar root find, with error propagation.
    abstract interface
        subroutine residual_fun(x, f, ierr)
            import :: dp
            real(dp), intent(in) :: x
            real(dp), intent(out) :: f
            integer, intent(out) :: ierr
        end subroutine residual_fun
    end interface

    !> LAPACK general linear solver (LU with partial pivoting).
    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dgesv
    end interface

contains

    !> Gauss-Legendre nodes and weights on an arbitrary interval.
    !!
    !! Nodes are obtained by Newton iteration on the Legendre polynomial
    !! `P_n`, using the standard Chebyshev-like initial guess.  Accuracy is at
    !! the level of the working precision for every order used here.
    !!
    !! @param[in]  n  Number of nodes (must be >= 1)
    !! @param[in]  a  Lower limit of integration
    !! @param[in]  b  Upper limit of integration
    !! @param[out] x  Nodes, ascending, size `n`
    !! @param[out] w  Weights, size `n`
    subroutine gauss_legendre(n, a, b, x, w)
        integer, intent(in) :: n
        real(dp), intent(in) :: a, b
        real(dp), intent(out) :: x(:), w(:)

        integer :: i, j, iter
        real(dp) :: z, z_old, p0, p1, p2, dp1, half, mid

        if (n < 1) return
        if (size(x) < n .or. size(w) < n) return

        half = 0.5_dp * (b - a)
        mid = 0.5_dp * (b + a)

        do i = 1, (n + 1) / 2
            z = cos(PI * (real(i, dp) - 0.25_dp) / (real(n, dp) + 0.5_dp))
            do iter = 1, 100
                ! Legendre recurrence: p2 = P_n(z), p1 = P_{n-1}(z).
                p0 = 1.0_dp
                p1 = 1.0_dp
                p2 = z
                do j = 2, n
                    p0 = p1
                    p1 = p2
                    p2 = (real(2 * j - 1, dp) * z * p1 - real(j - 1, dp) * p0) &
                         / real(j, dp)
                end do
                if (n == 1) then
                    dp1 = 1.0_dp
                else
                    dp1 = real(n, dp) * (z * p2 - p1) / (z * z - 1.0_dp)
                end if
                z_old = z
                z = z_old - p2 / dp1
                if (abs(z - z_old) <= 1.0e-16_dp * max(1.0_dp, abs(z))) exit
            end do

            x(i) = mid - half * z
            x(n + 1 - i) = mid + half * z
            w(i) = 2.0_dp * half / ((1.0_dp - z * z) * dp1 * dp1)
            w(n + 1 - i) = w(i)
        end do
    end subroutine gauss_legendre

    !> Dyadically graded Gauss-Legendre mesh on `[0, B]`.
    !!
    !! The spin kernels `a1` and `a2` have width `O(u)` while `B` ranges over
    !! many decades, so panels double in length away from the origin. The final
    !! non-dyadic interval is subdivided into pieces no wider than `1.5 u`.
    !! This keeps the kernel resolved everywhere with `O(log B)` panels while
    !! avoiding a U-dependent residual-panel geometry.
    !!
    !! @param[in]  B             Upper limit (must be > 0)
    !! @param[in]  u             Kernel scale `U/4`
    !! @param[in]  n_per_panel   Gauss-Legendre nodes per panel
    !! @param[out] x             Nodes (allocated here on success)
    !! @param[out] w             Weights (allocated here on success)
    !! @param[out] ierr          `ERROR_SUCCESS` or `ERROR_CONVERGENCE_FAILED`
    subroutine lambda_mesh(B, u, n_per_panel, x, w, ierr)
        real(dp), intent(in) :: B, u
        integer, intent(in) :: n_per_panel
        real(dp), allocatable, intent(out) :: x(:), w(:)
        integer, intent(out) :: ierr

        real(dp) :: edges(MAX_LAMBDA_PANELS + 1), v, scale, residual_start
        integer :: n_edges, n_residual, p, i0
        real(dp), allocatable :: xp(:), wp(:)

        ierr = ERROR_SUCCESS

        ! The first panel resolves the kernel width itself; no floor is applied
        ! because `U_QUAD_MIN` already keeps `u` away from zero, and a floor
        ! here would only hide the reason for that limit.
        scale = u
        edges(1) = 0.0_dp
        n_edges = 1
        v = min(0.125_dp * scale, 0.25_dp * B)
        do while (v < 0.75_dp * B .and. n_edges < MAX_LAMBDA_PANELS)
            n_edges = n_edges + 1
            edges(n_edges) = v
            v = 2.0_dp * v
        end do
        residual_start = edges(n_edges)
        n_residual = max(1, ceiling((B - residual_start) / &
                                    (RESIDUAL_PANEL_MAX_WIDTH * scale)))
        if (n_edges + n_residual > size(edges)) then
            ierr = ERROR_CONVERGENCE_FAILED
            return
        end if
        do p = 1, n_residual
            n_edges = n_edges + 1
            edges(n_edges) = residual_start + real(p, dp) * (B - residual_start) / real(n_residual, dp)
        end do

        allocate(x((n_edges - 1) * n_per_panel))
        allocate(w((n_edges - 1) * n_per_panel))
        allocate(xp(n_per_panel), wp(n_per_panel))

        do p = 1, n_edges - 1
            call gauss_legendre(n_per_panel, edges(p), edges(p + 1), xp, wp)
            i0 = (p - 1) * n_per_panel
            x(i0 + 1:i0 + n_per_panel) = xp
            w(i0 + 1:i0 + n_per_panel) = wp
        end do

        deallocate(xp, wp)
    end subroutine lambda_mesh

    !> Number of `Lambda` nodes per panel below which the spin sector is
    !! under-resolved at small magnetization, as a function of `u = U / 4`.
    !!
    !! The XC potentials are finite differences of `e_xc`, and next to `m = 0`
    !! the signal they have to resolve is `e(n, m) - e(n, 0) = O(K m^2)`, i.e.
    !! 1e-11 of `e_xc` at `m / n = 1e-5`.  Without a floor the default order
    !! returned pure quadrature noise there, which flipped the **sign** of
    !! `V_xc_up - V_xc_dn` as `m -> 0` - the quantity that decides the Stoner
    !! instability of the SCF.
    !!
    !! The requirement is **not** a monotone function of `U`, which is why this
    !! is an envelope and not a fit.  The graded `Lambda` mesh ends with a
    !! leftover panel `[v, B]` whose length relative to the kernel width `u`
    !! depends on where the dyadic doubling happens to land, i.e. on the
    !! fractional part of `log2(B / u)`.  The order needed therefore oscillates
    !! with `U` instead of decreasing with it.
    !!
    !! Measured relative error of the `m -> 0` exchange splitting against
    !! `n_lambda = 80`, worst case over `n` in {0.3, 0.8, 0.95} and `m/n` in
    !! {1e-5, 1e-4}, on `U` from 0.5 to 2.4 in steps of 0.025 to 0.05
    !! ("sign" = the splitting came out with the wrong sign):
    !!
    !! | `U`    | 12    | 16    | 20    | 24   | 28    | 32   | 36   | 40   |
    !! |--------|-------|-------|-------|------|-------|------|------|------|
    !! | 0.65   |       |       |       |      | sign  | 7.7% | 0.5% | 0.1% |
    !! | 0.85   |       |       |       | sign | 38%   | 0.9% | 0.0% | 0.0% |
    !! | 0.90   |       |       |       | 52%  | 2.3%  | 0.1% | 0.0% | 0.0% |
    !! | 1.05   |       | sign  | 6.3%  | 0.0% | 0.0%  | 0.0% | 0.0% | 0.0% |
    !! | 1.20   |       | 5.9%  | 0.0%  | 0.0% | 0.0%  | 0.0% | 0.0% | 0.0% |
    !! | 1.35   | 57%   | 0.0%  | 0.0%  | 0.0% | 0.0%  | 0.0% | 0.0% | 0.0% |
    !! | 1.50   | 2.1%  | 0.0%  | 0.0%  | 0.0% | 0.0%  | 0.0% | 0.0% | 0.0% |
    !! | 2.00   | 0.0%  | 0.0%  | 0.0%  | 0.0% | 0.0%  | 0.0% | 0.0% | 0.0% |
    !!
    !! Blank cells are orders that the *previous* floor already clamped upwards,
    !! so no independent measurement of them exists - and that clamping is what
    !! made the earlier calibration look better than it was.
    !!
    !! Each branch below carries the order required at the **hardest** `U` of
    !! its whole interval, not at one endpoint.  The previous ladder binned by
    !! the upper edge of each step while the requirement does not decrease
    !! monotonically, so `U = 1.02` fell into the step calibrated for `U = 1.5`
    !! and came out with the splitting sign inverted.  End to end, with the
    !! envelope below and the production quadrature, the worst relative error
    !! over the same `(n, m/n)` set on `U` from 0.5 to 2.4 in steps of 0.025 is
    !! 0.1%, with no sign inversion anywhere.
    !!
    !! Above `U = 2` the nominal default of 12 was measured sufficient up to
    !! `U = 8.0` (same `(n, m/n)` set, `U` from 2.4 to 8.0: worst relative
    !! error 0.0210%, no sign inversion). The reference refines all three
    !! quadrature dimensions (`n_lambda = 64`, `n_k = 96`, `n_omega = 64`), as
    !! does `test_low_m_quadrature_floor_convergence`. Above `U = 8` it is NOT measured -
    !! tables are shipped up to `U = 20` - and since the requirement oscillates
    !! with `frac(log2(B / u))` rather than decreasing with `U`, the measured
    !! range must not be extrapolated.  Over that measured range the strongly
    !! coupled tables keep their established cost; the price of the envelope is
    !! paid by the weak and intermediate coupling range, and only next to
    !! `m = 0` (see `n_lambda_per_panel`).
    !!
    !! @param[in] u  Kernel scale `U / 4`
    !! @return       Minimum number of Gauss-Legendre nodes per `Lambda` panel
    pure function n_lambda_floor(u) result(n_min)
        real(dp), intent(in) :: u
        integer :: n_min

        if (u < 0.25_dp) then
            ! U < 1: the oscillation peaks here; U = 0.65 and U = 0.85 invert
            ! the sign of the splitting at every order up to 28.
            n_min = 40
        else if (u < 0.5_dp) then
            ! 1 <= U < 2.  The drop to the nominal order must wait until here:
            ! at U = 1.35 the nominal 12 misses the splitting by 57%.
            n_min = 28
        else
            ! U >= 2: the nominal order is converged over the measured range
            ! (2.4 <= U <= 8.0, worst error 0.0210%); U > 8 is unmeasured.
            n_min = 12
        end if
    end function n_lambda_floor

    !> Per-panel `Lambda` order actually used at a given spin cut-off `B`.
    !!
    !! `n_lambda_floor` is only needed where the signal is the `O(K m^2)`
    !! curvature next to `m = 0`, and `m ~ n exp(-pi B / (2 u))` turns that into
    !! a condition on `B / u` alone (see `B_FLOOR_FRAC`, which carries the
    !! measured cost/accuracy curve).  At larger `m` the nominal order is
    !! already converged, and skipping the floor there is what keeps a weak
    !! coupling table affordable: the order enters the `n_k + n_Lambda` linear
    !! solve cubically.  Note that the `(Q, B)` inversion visits large `B` for
    !! *every* grid point, so the saving is 24% at `U = 1`, not a factor.
    !!
    !! @param[in] u          Kernel scale `U / 4`
    !! @param[in] B          Spin rapidity cut-off
    !! @param[in] n_nominal  Requested nodes per panel
    !! @return               Nodes per panel to use
    pure function n_lambda_per_panel(u, B, n_nominal) result(n_per)
        real(dp), intent(in) :: u, B
        integer, intent(in) :: n_nominal
        integer :: n_per

        n_per = max(4, n_nominal)
        if (B >= B_FLOOR_FRAC * u) n_per = max(n_per, n_lambda_floor(u))
    end function n_lambda_per_panel

    !> Largest useful spin cut-off for a given kernel scale `u = U/4`.
    !!
    !! The magnetization vanishes as `m ~ exp(-pi B / (2 u))`, so there is no
    !! point in searching beyond a `B` where `m` has already underflowed the
    !! working precision.  Keeping the bound tight also keeps the graded
    !! `Lambda` mesh short, which dominates the cost per linear solve.
    pure function b_max_for(u) result(b_max)
        real(dp), intent(in) :: u
        real(dp) :: b_max

        ! B/u = 25 makes the omitted magnetization exponentially smaller than
        ! double precision in the asymptotic estimate. A fixed B=30 floor was
        ! needlessly expensive at weak coupling (U=1 used B=30 although
        ! 25*u=6.25 already lies beyond the numerical information content).
        b_max = max(25.0_dp * u, 10.0_dp * B_MIN)
    end function b_max_for

    !> Initial guess for the spin cut-off from `m ~ (n) exp(-pi B / (2 u))`.
    pure function b_guess_for(n_tot, m, u) result(b0)
        real(dp), intent(in) :: n_tot, m, u
        real(dp) :: b0

        b0 = (2.0_dp * u / PI) * log(1.0_dp + n_tot / max(m, 1.0e-300_dp))
        b0 = min(max(b0, 10.0_dp * B_MIN), b_max_for(u))
    end function b_guess_for

    !> Lieb-Wu kernel `a_j(x) = (1/2 pi) 2 j u / ((j u)^2 + x^2)`.
    pure function a_kernel(j, u, x) result(a)
        integer, intent(in) :: j
        real(dp), intent(in) :: u, x
        real(dp) :: a, ju

        ju = real(j, dp) * u
        a = ju / (PI * (ju * ju + x * x))
    end function a_kernel

    !> Solve the coupled Lieb-Wu integral equations at fixed `(Q, B)`.
    !!
    !! The parity of `rho` and `sigma` is used to fold both integrals onto the
    !! positive half axis, so the linear system has size `n_k + n_Lambda`
    !! rather than twice that.
    !!
    !! @param[in]  Q      Charge rapidity cut-off, in `(0, pi]`
    !! @param[in]  B      Spin rapidity cut-off, `>= 0`
    !! @param[in]  U      Hubbard interaction, must be positive
    !! @param[in]  quad   Quadrature orders
    !! @param[out] n_tot  Total density `n_up + n_down`
    !! @param[out] n_dn   Minority (spin-down) density
    !! @param[out] e_ba   Bethe Ansatz total energy per site
    !! @param[out] ierr   `ERROR_SUCCESS`, or a LAPACK/validation error code
    subroutine lieb_wu_solve_qb(Q, B, U, quad, n_tot, n_dn, e_ba, ierr)
        real(dp), intent(in) :: Q, B, U
        type(lw_quad_t), intent(in) :: quad
        real(dp), intent(out) :: n_tot, n_dn, e_ba
        integer, intent(out) :: ierr

        real(dp), allocatable :: xk(:), wk(:), sk(:), ck(:)
        real(dp), allocatable :: xl(:), wl(:)
        real(dp), allocatable :: A(:, :), rhs(:, :)
        integer, allocatable :: ipiv(:)
        real(dp) :: u4
        integer :: nk, nl, ntot_eq, i, j, info

        ierr = ERROR_SUCCESS
        n_tot = 0.0_dp
        n_dn = 0.0_dp
        e_ba = 0.0_dp

        if (U <= 0.0_dp .or. Q < 0.0_dp .or. Q > PI + N_FULL_TOL .or. B < 0.0_dp) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (Q <= 0.0_dp) return

        u4 = 0.25_dp * U

        ! Fully polarized limit: sigma vanishes and rho is flat.
        if (B <= B_MIN) then
            n_tot = Q / PI
            n_dn = 0.0_dp
            e_ba = -2.0_dp * sin(Q) / PI
            return
        end if

        nk = max(4, quad%n_k)
        allocate(xk(nk), wk(nk), sk(nk), ck(nk))
        call gauss_legendre(nk, 0.0_dp, Q, xk, wk)
        sk = sin(xk)
        ck = cos(xk)

        call lambda_mesh(B, u4, n_lambda_per_panel(u4, B, quad%n_lambda), xl, wl, ierr)
        if (ierr /= ERROR_SUCCESS) return
        nl = size(xl)
        ntot_eq = nk + nl

        allocate(A(ntot_eq, ntot_eq), rhs(ntot_eq, 1), ipiv(ntot_eq))
        A = 0.0_dp
        rhs = 0.0_dp

        ! Charge equation, rows 1..nk.
        do i = 1, nk
            A(i, i) = 1.0_dp
            do j = 1, nl
                A(i, nk + j) = -ck(i) * wl(j) &
                    * (a_kernel(1, u4, sk(i) - xl(j)) + a_kernel(1, u4, sk(i) + xl(j)))
            end do
            rhs(i, 1) = 1.0_dp / TWOPI
        end do

        ! Spin equation, rows nk+1..nk+nl.
        do i = 1, nl
            do j = 1, nk
                A(nk + i, j) = -wk(j) &
                    * (a_kernel(1, u4, xl(i) - sk(j)) + a_kernel(1, u4, xl(i) + sk(j)))
            end do
            do j = 1, nl
                A(nk + i, nk + j) = A(nk + i, nk + j) + wl(j) &
                    * (a_kernel(2, u4, xl(i) - xl(j)) + a_kernel(2, u4, xl(i) + xl(j)))
            end do
            A(nk + i, nk + i) = A(nk + i, nk + i) + 1.0_dp
        end do

        call dgesv(ntot_eq, 1, A, ntot_eq, ipiv, rhs, ntot_eq, info)
        if (info > 0) then
            ierr = ERROR_SINGULAR_MATRIX
        else if (info < 0) then
            ierr = ERROR_LAPACK_INVALID_ARG
        else
            n_tot = 2.0_dp * dot_product(wk, rhs(1:nk, 1))
            n_dn = 2.0_dp * dot_product(wl, rhs(nk + 1:ntot_eq, 1))
            e_ba = -4.0_dp * dot_product(wk * ck, rhs(1:nk, 1))
            if (.not. (ieee_is_finite(n_tot) .and. ieee_is_finite(n_dn) &
                       .and. ieee_is_finite(e_ba))) ierr = ERROR_NOT_A_NUMBER
        end if

        deallocate(xk, wk, sk, ck, xl, wl, A, rhs, ipiv)
    end subroutine lieb_wu_solve_qb

    !> Solve the Lieb-Wu equations in the unpolarized limit `B -> infinity`.
    !!
    !! Eliminating `sigma` by Fourier transform on the whole line leaves
    !!
    !! ```
    !! rho(k) = 1/(2 pi) + cos(k) Int_{-Q}^{Q} R(sin k - sin k') rho(k') dk'
    !! R(x)   = (1/pi) Int_0^inf cos(w x) / (1 + exp(2 u w)) dw
    !! ```
    !!
    !! Folding onto `[0, Q]` gives `R(s - s') + R(s + s') = (2/pi) Int_0^inf
    !! cos(w s) cos(w s') / (1 + exp(2 u w)) dw`, a single separable Fourier
    !! quadrature that is assembled as one matrix product.
    !!
    !! @param[in]  Q      Charge rapidity cut-off, in `(0, pi]`
    !! @param[in]  U      Hubbard interaction, must be positive
    !! @param[in]  quad   Quadrature orders
    !! @param[out] n_tot  Total density (the minority density is `n_tot/2`)
    !! @param[out] e_ba   Bethe Ansatz total energy per site
    !! @param[out] ierr   `ERROR_SUCCESS`, or a LAPACK/validation error code
    subroutine lieb_wu_solve_unpolarized(Q, U, quad, n_tot, e_ba, ierr)
        real(dp), intent(in) :: Q, U
        type(lw_quad_t), intent(in) :: quad
        real(dp), intent(out) :: n_tot, e_ba
        integer, intent(out) :: ierr

        real(dp), allocatable :: xk(:), wk(:), sk(:), ck(:)
        real(dp), allocatable :: xw(:), ww(:)
        real(dp), allocatable :: cmat(:, :), cmat_d(:, :), ksym(:, :)
        real(dp), allocatable :: A(:, :), rhs(:, :)
        integer, allocatable :: ipiv(:)
        real(dp) :: u4, omega_max, panel
        integer :: nk, nw, n_panels, n_per, i, j, p, i0, info

        ierr = ERROR_SUCCESS
        n_tot = 0.0_dp
        e_ba = 0.0_dp

        if (U <= 0.0_dp .or. Q < 0.0_dp .or. Q > PI + N_FULL_TOL) then
            ierr = ERROR_INVALID_INPUT
            return
        end if
        if (Q <= 0.0_dp) return

        u4 = 0.25_dp * U

        nk = max(4, quad%n_k)
        allocate(xk(nk), wk(nk), sk(nk), ck(nk))
        call gauss_legendre(nk, 0.0_dp, Q, xk, wk)
        sk = sin(xk)
        ck = cos(xk)

        ! omega mesh: uniform panels up to the exponential cut-off.  Panel
        ! length is kept near one so that cos(omega sin k) stays resolved.
        omega_max = OMEGA_CUT / (2.0_dp * u4)
        n_panels = min(MAX_OMEGA_PANELS, max(8, ceiling(omega_max)))
        n_per = max(4, quad%n_omega)
        nw = n_panels * n_per
        panel = omega_max / real(n_panels, dp)

        allocate(xw(nw), ww(nw))
        do p = 1, n_panels
            i0 = (p - 1) * n_per
            call gauss_legendre(n_per, real(p - 1, dp) * panel, real(p, dp) * panel, &
                                xw(i0 + 1:i0 + n_per), ww(i0 + 1:i0 + n_per))
        end do

        allocate(cmat(nk, nw), cmat_d(nk, nw), ksym(nk, nk))
        do p = 1, nw
            do i = 1, nk
                cmat(i, p) = cos(xw(p) * sk(i))
            end do
            cmat_d(:, p) = cmat(:, p) * ww(p) / (1.0_dp + exp(2.0_dp * u4 * xw(p)))
        end do
        ksym = (2.0_dp / PI) * matmul(cmat_d, transpose(cmat))

        allocate(A(nk, nk), rhs(nk, 1), ipiv(nk))
        do i = 1, nk
            do j = 1, nk
                A(i, j) = -ck(i) * ksym(i, j) * wk(j)
            end do
            A(i, i) = A(i, i) + 1.0_dp
            rhs(i, 1) = 1.0_dp / TWOPI
        end do

        call dgesv(nk, 1, A, nk, ipiv, rhs, nk, info)
        if (info > 0) then
            ierr = ERROR_SINGULAR_MATRIX
        else if (info < 0) then
            ierr = ERROR_LAPACK_INVALID_ARG
        else
            n_tot = 2.0_dp * dot_product(wk, rhs(:, 1))
            e_ba = -4.0_dp * dot_product(wk * ck, rhs(:, 1))
            if (.not. (ieee_is_finite(n_tot) .and. ieee_is_finite(e_ba))) &
                ierr = ERROR_NOT_A_NUMBER
        end if

        deallocate(xk, wk, sk, ck, xw, ww, cmat, cmat_d, ksym, A, rhs, ipiv)
    end subroutine lieb_wu_solve_unpolarized

    !> Bracketed root find with a warm start (Illinois / modified false position).
    !!
    !! A bracket is grown geometrically around `x_guess` until the residual
    !! changes sign, then Illinois iterations narrow it.  The bracket is never
    !! left, so the iteration cannot diverge; a wrong warm start only costs a
    !! few extra residual evaluations.
    !!
    !! @param[in]  fun       Residual, must be continuous and monotone on `[x_lo, x_hi]`
    !! @param[in]  x_lo      Lower bound of the search interval
    !! @param[in]  x_hi      Upper bound of the search interval
    !! @param[in]  x_guess   Warm start (clamped into the bounds)
    !! @param[in]  span0     Initial half-width of the bracket around `x_guess`
    !! @param[in]  tol       Residual tolerance
    !! @param[in]  max_iter  Maximum Illinois iterations
    !! @param[out] x         Root
    !! @param[out] ierr      `ERROR_SUCCESS`, `ERROR_CONVERGENCE_FAILED`, or the
    !!                       error code raised by `fun`
    recursive subroutine solve_root_1d(fun, x_lo, x_hi, x_guess, span0, tol, &
                                       max_iter, x, ierr)
        procedure(residual_fun) :: fun
        real(dp), intent(in) :: x_lo, x_hi, x_guess, span0, tol
        integer, intent(in) :: max_iter
        real(dp), intent(out) :: x
        integer, intent(out) :: ierr

        real(dp) :: a, b, fa, fb, fx, span, a_new, b_new
        integer :: it, side

        x = min(max(x_guess, x_lo), x_hi)
        span = max(span0, 1.0e-12_dp)

        a = max(x_lo, x - span)
        b = min(x_hi, x + span)
        if (b <= a) then
            a = x_lo
            b = x_hi
        end if

        call fun(a, fa, ierr)
        if (ierr /= ERROR_SUCCESS) return
        call fun(b, fb, ierr)
        if (ierr /= ERROR_SUCCESS) return

        do it = 1, 60
            if (fa * fb <= 0.0_dp) exit
            if (a <= x_lo .and. b >= x_hi) exit
            span = 3.0_dp * span
            a_new = max(x_lo, x - span)
            b_new = min(x_hi, x + span)
            if (a_new < a) then
                a = a_new
                call fun(a, fa, ierr)
                if (ierr /= ERROR_SUCCESS) return
            end if
            if (b_new > b) then
                b = b_new
                call fun(b, fb, ierr)
                if (ierr /= ERROR_SUCCESS) return
            end if
        end do

        if (fa * fb > 0.0_dp) then
            ! No sign change anywhere in the admissible interval.  Accept an
            ! endpoint only if it already satisfies the tolerance.
            if (abs(fa) <= tol) then
                x = a
                call fun(x, fx, ierr)
            else if (abs(fb) <= tol) then
                x = b
                call fun(x, fx, ierr)
            else
                ierr = ERROR_CONVERGENCE_FAILED
            end if
            return
        end if

        if (abs(fa) <= tol) then
            x = a
            call fun(x, fx, ierr)
            return
        end if
        if (abs(fb) <= tol) then
            x = b
            call fun(x, fx, ierr)
            return
        end if

        side = 0
        do it = 1, max_iter
            if (abs(fb - fa) <= tiny(1.0_dp)) exit
            x = (fb * a - fa * b) / (fb - fa)
            x = min(max(x, a), b)
            call fun(x, fx, ierr)
            if (ierr /= ERROR_SUCCESS) return
            if (abs(fx) <= tol) return
            if (b - a <= 1.0e-15_dp * max(1.0_dp, abs(x))) return
            if (fx * fb > 0.0_dp) then
                b = x
                fb = fx
                if (side == -1) fa = 0.5_dp * fa
                side = -1
            else
                a = x
                fa = fx
                if (side == 1) fb = 0.5_dp * fb
                side = 1
            end if
        end do

        ierr = ERROR_CONVERGENCE_FAILED
    end subroutine solve_root_1d

    !> Invert `(n, n_down) -> (Q, B)` and return the Bethe Ansatz energy.
    !!
    !! The inversion is nested: the outer root find fixes `w = atan(B)` from
    !! the minority density, and for every trial `w` an inner root find fixes
    !! `Q` from the total density.  Both residuals are monotone, so both levels
    !! are bracketed and cannot fail to converge on physical input.
    !!
    !! Special branches, in order of priority: empty band, full polarization
    !! (`n_down = 0`, `B = 0`), vanishing magnetization (`B = infinity`, Shiba
    !! reduction), and half filling (`Q = pi` exactly, no inversion in `Q`).
    !!
    !! @param[in]     n_tot   Total density, in `(0, 1]`
    !! @param[in]     n_dn_t  Minority density, in `[0, n_tot/2]`
    !! @param[in]     U       Hubbard interaction, must be positive
    !! @param[in]     quad    Quadrature orders and tolerances
    !! @param[out]    Q       Converged charge rapidity cut-off
    !! @param[out]    w       Converged spin angle, `B = tan(w)`, `pi/2` if `m = 0`
    !! @param[out]    e_ba    Bethe Ansatz total energy per site
    !! @param[out]    ierr    `ERROR_SUCCESS` or an error code
    !! @param[inout]  seed    Optional warm start; updated on success
    subroutine lieb_wu_invert(n_tot, n_dn_t, U, quad, Q, w, e_ba, ierr, seed)
        real(dp), intent(in) :: n_tot, n_dn_t, U
        type(lw_quad_t), intent(in) :: quad
        real(dp), intent(out) :: Q, w, e_ba
        integer, intent(out) :: ierr
        type(lw_seed_t), intent(inout), optional :: seed

        real(dp) :: m, B_cur, Q_cur, n_cur, n_dn_cur, e_cur
        real(dp) :: Q_guess, w_guess, dummy, b_max
        logical :: fix_Q

        ierr = ERROR_SUCCESS
        Q = 0.0_dp
        w = 0.0_dp
        e_ba = 0.0_dp
        B_cur = 0.0_dp
        Q_cur = 0.0_dp
        n_cur = 0.0_dp
        n_dn_cur = 0.0_dp
        e_cur = 0.0_dp

        if (U < U_QUAD_MIN) then
            ierr = ERROR_INVALID_INPUT
            return
        end if
        if (n_tot < -N_MINORITY_TOL .or. n_tot > 1.0_dp + N_FULL_TOL &
            .or. n_dn_t < -N_MINORITY_TOL &
            .or. n_dn_t > 0.5_dp * n_tot + N_FULL_TOL) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (n_tot <= N_MINORITY_TOL) return

        m = n_tot - 2.0_dp * n_dn_t

        ! Fully polarized: no spin rapidities, free spinless band.
        if (n_dn_t <= N_MINORITY_TOL) then
            Q = PI * n_tot
            w = 0.0_dp
            e_ba = -2.0_dp * sin(Q) / PI
            if (present(seed)) then
                seed%valid = .true.
                seed%Q = Q
                seed%w = w
            end if
            return
        end if

        ! Unpolarized: B = infinity, Shiba reduction, single unknown Q.
        if (m <= M_ZERO_TOL) then
            w = 0.5_dp * PI
            if (n_tot >= 1.0_dp - N_FULL_TOL) then
                Q = PI
                call res_q_unpol(Q, dummy, ierr)
            else
                Q_guess = PI * n_tot
                if (present(seed)) then
                    if (seed%valid) Q_guess = seed%Q
                end if
                call solve_root_1d(res_q_unpol, 1.0e-10_dp, PI, Q_guess, &
                                   0.05_dp * PI, quad%tol, quad%max_iter, Q, ierr)
                if (ierr == ERROR_SUCCESS) call res_q_unpol(Q, dummy, ierr)
            end if
            if (ierr /= ERROR_SUCCESS) return
            e_ba = e_cur
            if (present(seed)) then
                seed%valid = .true.
                seed%Q = Q
                seed%w = w
            end if
            return
        end if

        ! General case.  Half filling pins Q = pi: integrating the charge
        ! equation over the full band gives n = 1 for any B.
        fix_Q = (n_tot >= 1.0_dp - N_FULL_TOL)

        b_max = b_max_for(0.25_dp * U)
        Q_guess = 0.25_dp * PI * (3.0_dp * n_tot + m)
        w_guess = atan(b_guess_for(n_tot, m, 0.25_dp * U))
        if (present(seed)) then
            if (seed%valid) then
                Q_guess = seed%Q
                if (seed%w > 0.0_dp .and. seed%w < 0.5_dp * PI) w_guess = seed%w
            end if
        end if
        Q_cur = min(max(Q_guess, 1.0e-10_dp), PI)

        call solve_root_1d(res_w, atan(B_MIN), atan(b_max), w_guess, &
                           0.05_dp, quad%tol, quad%max_iter, w, ierr)
        if (ierr /= ERROR_SUCCESS) return

        ! Refresh the cached inner solution so that it belongs to the root.
        call res_w(w, dummy, ierr)
        if (ierr /= ERROR_SUCCESS) return

        Q = Q_cur
        e_ba = e_cur
        if (present(seed)) then
            seed%valid = .true.
            seed%Q = Q
            seed%w = w
        end if

    contains

        !> Total-density residual of the unpolarized (B = infinity) branch.
        subroutine res_q_unpol(x, f, jerr)
            real(dp), intent(in) :: x
            real(dp), intent(out) :: f
            integer, intent(out) :: jerr

            call lieb_wu_solve_unpolarized(x, U, quad, n_cur, e_cur, jerr)
            n_dn_cur = 0.5_dp * n_cur
            f = n_cur - n_tot
        end subroutine res_q_unpol

        !> Total-density residual at the cached spin cut-off `B_cur`.
        subroutine res_q(x, f, jerr)
            real(dp), intent(in) :: x
            real(dp), intent(out) :: f
            integer, intent(out) :: jerr

            call lieb_wu_solve_qb(x, B_cur, U, quad, n_cur, n_dn_cur, e_cur, jerr)
            f = n_cur - n_tot
        end subroutine res_q

        !> Minority-density residual; each evaluation re-solves for `Q`.
        subroutine res_w(x, f, jerr)
            real(dp), intent(in) :: x
            real(dp), intent(out) :: f
            integer, intent(out) :: jerr

            real(dp) :: q_root, ignored

            B_cur = tan(x)
            if (fix_Q) then
                Q_cur = PI
                call res_q(PI, ignored, jerr)
            else
                call solve_root_1d(res_q, 1.0e-10_dp, PI, Q_cur, 0.02_dp * PI, &
                                   quad%tol, quad%max_iter, q_root, jerr)
                if (jerr == ERROR_SUCCESS) then
                    Q_cur = q_root
                    call res_q(Q_cur, ignored, jerr)
                end if
            end if
            f = n_dn_cur - n_dn_t
        end subroutine res_w

    end subroutine lieb_wu_invert

    !> Bethe Ansatz total energy per site at density `n` and magnetization `m`.
    !!
    !! Thin wrapper over `lieb_wu_invert`, kept for validation against closed
    !! forms (half filling, strong coupling, free limit).
    !!
    !! @param[in]  n     Total density, in `(0, 1]`
    !! @param[in]  m     Magnetization `n_up - n_down`, in `[0, n]`
    !! @param[in]  U     Hubbard interaction, must be positive
    !! @param[in]  quad  Quadrature orders and tolerances
    !! @param[out] e_ba  Total energy per site
    !! @param[out] ierr  `ERROR_SUCCESS` or an error code
    subroutine lieb_wu_energy(n, m, U, quad, e_ba, ierr)
        real(dp), intent(in) :: n, m, U
        type(lw_quad_t), intent(in) :: quad
        real(dp), intent(out) :: e_ba
        integer, intent(out) :: ierr

        real(dp) :: Q, w

        call lieb_wu_invert(n, 0.5_dp * (n - abs(m)), U, quad, Q, w, e_ba, ierr)
    end subroutine lieb_wu_energy

    !> Exchange-correlation energy per site, `e_xc = e_BA - T_s - U n_up n_dn`.
    !!
    !! Spin exchange (`m < 0`) and particle-hole symmetry (`n > 1`) are applied
    !! first, so any admissible `(n_up, n_dn)` in `[0, 1]^2` is accepted.  The
    !! result is the quantity stored in the `exc` column of the reference
    !! tables: the Hartree term `U n_up n_dn` and the Kohn-Sham kinetic energy
    !! `T_s` of the thermodynamic limit are both already subtracted.
    !!
    !! @param[in]     n_up  Spin-up density, in `[0, 1]`
    !! @param[in]     n_dn  Spin-down density, in `[0, 1]`
    !! @param[in]     U     Hubbard interaction, must be positive
    !! @param[in]     quad  Quadrature orders and tolerances
    !! @param[out]    exc   XC energy per site
    !! @param[out]    ierr  `ERROR_SUCCESS` or an error code
    !! @param[inout]  seed  Optional warm start for the `(Q, B)` inversion
    subroutine lieb_wu_exc(n_up, n_dn, U, quad, exc, ierr, seed)
        real(dp), intent(in) :: n_up, n_dn, U
        type(lw_quad_t), intent(in) :: quad
        real(dp), intent(out) :: exc
        integer, intent(out) :: ierr
        type(lw_seed_t), intent(inout), optional :: seed

        real(dp) :: n_maj, n_min, swap, n_total, Q, w, e_ba, t_s

        exc = 0.0_dp
        ierr = ERROR_SUCCESS

        if (abs(U) < U_SMALL) return
        if (U < U_QUAD_MIN) then
            ierr = ERROR_INVALID_INPUT
            return
        end if
        if (n_up < -N_MINORITY_TOL .or. n_dn < -N_MINORITY_TOL &
            .or. n_up > 1.0_dp + N_INPUT_TOL .or. n_dn > 1.0_dp + N_INPUT_TOL) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        n_maj = max(n_up, n_dn)
        n_min = min(n_up, n_dn)

        ! Particle-hole symmetry: e_xc(n_up, n_dn) = e_xc(1-n_up, 1-n_dn).
        if (n_maj + n_min > 1.0_dp) then
            swap = 1.0_dp - n_min
            n_min = 1.0_dp - n_maj
            n_maj = swap
        end if

        n_maj = min(max(n_maj, 0.0_dp), 1.0_dp)
        n_min = min(max(n_min, 0.0_dp), n_maj)
        n_total = n_maj + n_min

        if (n_total <= N_MINORITY_TOL) return
        if (n_min <= N_MINORITY_TOL) return

        call lieb_wu_invert(n_total, n_min, U, quad, Q, w, e_ba, ierr, seed)
        if (ierr /= ERROR_SUCCESS) return

        t_s = -(2.0_dp / PI) * (sin(PI * n_maj) + sin(PI * n_min))
        exc = e_ba - t_s - U * n_maj * n_min
    end subroutine lieb_wu_exc

end module lieb_wu_integral
