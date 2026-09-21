!> Module for LSDA exchange-correlation functional evaluation
!!
!! This module provides a high-level interface for evaluating the XC functional
!! of the 1D Hubbard model at arbitrary points (n_up, n_dw) using:
!! 1. Pre-computed Bethe Ansatz tables
!! 2. 2D bicubic spline interpolation
!! 3. Physical symmetries to cover the full domain
!!
!! Physical symmetries (4 regions):
!!   Region I   (m ≥ 0, n ≤ 1): Direct table lookup
!!   Region II  (m < 0, n ≤ 1): Spin exchange symmetry
!!   Region III (m < 0, n > 1): Particle-hole symmetry
!!   Region IV  (m ≥ 0, n > 1): Combined symmetry
!!
!! Attractive interaction (U < 0) is handled by the Shiba (partial
!! particle-hole) transformation, applied on top of the four regions above:
!!   e_xc(n_up, n_dw; U<0)    =  e_xc(1-n_up, n_dw; |U|)
!!   V_xc^up(n_up, n_dw; U<0) = -V_xc^up(1-n_up, n_dw; |U|)
!!   V_xc^dn(n_up, n_dw; U<0) = +V_xc^dn(1-n_up, n_dw; |U|)
!! Only the |U| table is ever tabulated; the whole sign dependence enters
!! through these three relations, exactly as in the C++ reference
!! (`original/spline2D.cc`, functions exc_value / Vxc_up_value / Vxc_dn_value).
module xc_lsda
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use lsda_constants, only: dp, U_SMALL, PI, HALF_FILLING_SNAP_TOL
    use spline2d, only: spline2d_t, spline2d_init, spline2d_eval, spline2d_destroy
    use table_io, only: xc_table_t, read_fortran_table, deallocate_table
    use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT, ERROR_OUT_OF_BOUNDS, ERROR_FILE_READ, ERROR_SPLINE_INITIALIZATION_FAILED
    implicit none
    private

    !> LSDA XC functional type containing splines for exc, vxc_up, vxc_down
    type, public :: xc_lsda_t
        real(dp) :: U = 0.0_dp                  !< Hubbard U parameter, WITH sign (negative = attractive)
        real(dp) :: smoothing_width = 0.0_dp    !< Half-width w of the linear smoothing of V_xc around n = 1 (0 = off)
        type(spline2d_t) :: spl_exc             !< Spline for e_xc(n, m)
        type(spline2d_t) :: spl_vxc_up          !< Spline for V_xc^up(n, m)
        type(spline2d_t) :: spl_vxc_down        !< Spline for V_xc^dn(n, m)
        logical :: initialized = .false.        !< Initialization flag
    end type xc_lsda_t

    !> Largest admissible smoothing half-width (the window must stay inside 0 < n < 2)
    real(dp), parameter, public :: XC_SMOOTHING_WIDTH_MAX = 1.0_dp

    !> Tolerance used to match the signed interaction requested by a calculation
    !! against the interaction carried by the initialized XC functional.
    real(dp), parameter, public :: XC_U_MATCH_TOL = 1.0e-6_dp

    !> Width of the band in which a spin channel counts as empty
    !!
    !! Same constant as the C++ reference, where the tests read
    !! `fabs(dens - mag) < 1.0e-14` (empty spin-down) and
    !! `fabs(dens + mag) < 1.0e-14` (empty spin-up).
    real(dp), parameter :: EMPTY_CHANNEL_TOL = 1.0e-14_dp

    !> Width of the fully polarized corner shortcut n = m = 1 of the C++
    !! reference (`(1 - 1e-14 <= dens <= 1) && (1 - 1e-14 <= mag <= 1)`)
    real(dp), parameter :: CORNER_TOL = 1.0e-14_dp

    !> Half-width of the band in which (n, m) is snapped onto a region boundary
    !!
    !! The two region boundaries, n = n_up + n_dw = 1 (particle-hole line) and
    !! m = n_up - n_dw = 0 (spin-symmetry line), are the loci where the four
    !! symmetry branches meet, and the branches are NOT small perturbations of
    !! each other: crossing n = 1 flips the overall sign of V_xc, crossing m = 0
    !! exchanges the two spin channels. Roundoff of a few 1e-16 in n_up or n_dw
    !! (see Bug #1 of CLAUDE.md: at exact half filling n_up + n_dw may come out
    !! as 1.0000000000000002) must therefore not be allowed to decide the
    !! branch, which is why a tolerance is needed at all.
    real(dp), parameter :: REGION_SNAP_TOL = HALF_FILLING_SNAP_TOL

    !> Number of public XC evaluations since the last diagnostic reset.
    !!
    !! This is intentionally module-global diagnostic state, not physical
    !! state in `xc_lsda_t`: it lets regression tests prove that the SCF uses
    !! its output-density cache without changing the functional API.
    integer, save :: xc_evaluation_count = 0

    public :: xc_lsda_init
    public :: get_exc
    public :: get_vxc
    public :: xc_lsda_destroy
    public :: dexc_dndown_b0
    public :: reset_xc_evaluation_count, get_xc_evaluation_count

    private :: integral_1
    private :: is_corner_shortcut
    private :: exc_recursive
    private :: vxc_up_recursive
    private :: vxc_dn_recursive
    private :: determine_region
    private :: snap_to_boundaries
    private :: convert_to_nm
    private :: apply_symmetry_transform
    private :: eval_vxc_branch
    private :: get_exc_positive
    private :: get_vxc_positive

contains

    !> Reset the diagnostic count of public XC evaluations.
    subroutine reset_xc_evaluation_count()
        xc_evaluation_count = 0
    end subroutine reset_xc_evaluation_count

    !> Return the diagnostic count of public XC evaluations.
    function get_xc_evaluation_count() result(count)
        integer :: count
        count = xc_evaluation_count
    end function get_xc_evaluation_count

    !> Initialize XC functional from table file
    !!
    !! Loads table and constructs 2D splines for exc, vxc_up, vxc_down.
    !!
    !! The optional `smoothing_width` activates the linear smoothing of the
    !! V_xc discontinuity at n = 1 (see `get_vxc`). The default, 0, disables it
    !! and reproduces the C++ reference exactly.
    !!
    !! @param[out] xc         XC functional object
    !! @param[in]  table_file Optional path to table file (Fortran binary
    !!                         format). It is not required when `u_signed` is
    !!                         present and |U| < U_SMALL: the exact
    !!                         non-interacting functional has e_xc = V_xc = 0.
    !! @param[out] ierr Error code (0 = success)
    !! @param[in]  smoothing_width Optional half-width w of the V_xc smoothing
    !!                             window around n = 1; must be a non-NaN value
    !!                             with 0 <= w < XC_SMOOTHING_WIDTH_MAX
    !!                             (default 0 = off). NaN is rejected with
    !!                             ERROR_INVALID_INPUT because it would otherwise
    !!                             disable the smoothing silently in `get_vxc`.
    !! @param[in]  u_signed Optional Hubbard U of the run, WITH sign. Tables are
    !!                      only ever tabulated for |U|, so the table file alone
    !!                      cannot tell an attractive run from a repulsive one.
    !!                      Passing a negative value here activates the Shiba
    !!                      transformation in `get_exc`/`get_vxc`. The magnitude
    !!                      must agree with the table's |U| (otherwise the wrong
    !!                      table was loaded) or ERROR_INVALID_INPUT is returned.
    !!                      When absent, the sign from the table (positive) is
    !!                      kept and no Shiba transformation is applied.
    subroutine xc_lsda_init(xc, table_file, ierr, smoothing_width, u_signed)
        type(xc_lsda_t), intent(out) :: xc
        character(len=*), intent(in), optional :: table_file
        integer, intent(out) :: ierr
        real(dp), intent(in), optional :: smoothing_width
        real(dp), intent(in), optional :: u_signed

        type(xc_table_t) :: table
        integer :: io_stat
        integer, allocatable :: n_y_pts(:)
        real(dp), allocatable :: dexc_dm_first(:), dexc_dm_last(:)
        integer :: i

        ierr = ERROR_SUCCESS

        if (present(smoothing_width)) then
            ! NaN passes both range comparisons, is stored, and then fails
            ! SILENTLY: get_vxc evaluates `w > 0.0_dp` as false and takes the
            ! unsmoothed branch without any warning. Both infinities are already
            ! rejected by the range checks (+Inf >= 1 and -Inf < 0), so only NaN
            ! needs an explicit test here.
            if (ieee_is_nan(smoothing_width) .or. &
                smoothing_width < 0.0_dp .or. smoothing_width >= XC_SMOOTHING_WIDTH_MAX) then
                print *, "ERROR: xc smoothing width must be in [0, 1), got: ", smoothing_width
                ierr = ERROR_INVALID_INPUT
                return
            end if
            xc%smoothing_width = smoothing_width
        else
            xc%smoothing_width = 0.0_dp
        end if

        ! U = 0 is exactly the non-interacting gas.  Do not require a table:
        ! there is no exchange-correlation energy or potential to interpolate.
        ! This physical-zero decision is deliberately independent of
        ! XC_U_MATCH_TOL, which is only for matching a requested interaction
        ! to the interaction carried by an XC table.
        if (present(u_signed)) then
            if (abs(u_signed) < U_SMALL) then
                xc%U = 0.0_dp
                xc%initialized = .true.
                return
            end if
        end if

        if (.not. present(table_file)) then
            print *, "ERROR: an XC table file is required for nonzero U"
            ierr = ERROR_INVALID_INPUT
            return
        end if

        call read_fortran_table(table_file, table, io_stat)
        if (io_stat /= 0) then
            print *, "ERROR: Failed to read table file: ", trim(table_file)
            ierr = ERROR_FILE_READ
            return
        end if

        xc%U = table%U

        if (present(u_signed)) then
            ! The table is the one for |U|; the run's sign has to come from the
            ! caller. Reject a magnitude mismatch instead of silently running
            ! the functional of a different interaction strength.
            if (abs(abs(u_signed) - abs(table%U)) > XC_U_MATCH_TOL) then
                print *, "ERROR: xc table |U| does not match requested |U|: ", abs(table%U), abs(u_signed)
                ierr = ERROR_INVALID_INPUT
                call deallocate_table(table)
                return
            end if
            xc%U = sign(abs(table%U), u_signed)
        end if

        ! Create n_y_pts array (same for all n_grid points in current implementation)
        allocate(n_y_pts(table%n_points_n))
        do i = 1, table%n_points_n
            n_y_pts(i) = table%n_points_m
        end do

        ! Boundary conditions of the magnetization (row) splines, as in the C++
        ! reference (`original/spline2D.cc`, `build_spline_data`): clamped, not
        ! natural. A natural spline forces f'' = 0 at m = 0 and m = n, which is
        ! wrong at both ends of every row and shows up as an error of up to
        ! 3e-4 in e_xc close to full polarization.
        !
        ! For e_xc the two end slopes are known analytically:
        !  - at m = 0 spin symmetry gives e_xc(n, m) = e_xc(n, -m), hence
        !    ∂e_xc/∂m = 0;
        !  - at m = n the endpoint derivative follows from V_xc^dn(n, m = n),
        !    the last entry of the tabulated V_xc^dn row.
        !
        ! NOTE on the second one: along a fixed-n row, dn_up/dm = +1/2 and
        ! dn_down/dm = -1/2. Since e_xc vanishes identically on the fully
        ! polarized line, this gives ∂e_xc/∂m = -(1/2)V_xc^dn at m = n.
        ! original/spline2D.cc:378 passes -V_xc^dn despite stating the factor
        ! 1/2 in its own comment at lines 369-371. That is a known C++ bug: it
        ! doubles the endpoint slope and makes the last interval overshoot. We
        ! deliberately use the analytic derivative here rather than reproducing
        ! that error; parity with the reference is therefore not expected in the
        ! final magnetization interval closest to full polarization.
        ! For the two V_xc tables no closed form is available and the reference
        ! clamps with the secants of the first and last interval, which is the
        ! default of `spline2d_init` for `row_bc = 'clamped'`.
        allocate(dexc_dm_first(table%n_points_n), dexc_dm_last(table%n_points_n))
        do i = 1, table%n_points_n
            dexc_dm_first(i) = 0.0_dp
            dexc_dm_last(i) = -0.5_dp * table%vxc_down(n_y_pts(i), i)
        end do

        call spline2d_init(xc%spl_exc, table%n_grid, table%m_grid, &
                           table%exc, n_y_pts, row_bc = 'clamped', &
                           dy_first = dexc_dm_first, dy_last = dexc_dm_last)

        call spline2d_init(xc%spl_vxc_up, table%n_grid, table%m_grid, &
                           table%vxc_up, n_y_pts, row_bc = 'clamped')

        call spline2d_init(xc%spl_vxc_down, table%n_grid, table%m_grid, &
                           table%vxc_down, n_y_pts, row_bc = 'clamped')

        deallocate(n_y_pts, dexc_dm_first, dexc_dm_last)

        if (.not. xc%spl_exc%initialized .or. &
            .not. xc%spl_vxc_up%initialized .or. &
            .not. xc%spl_vxc_down%initialized) then
            ierr = ERROR_SPLINE_INITIALIZATION_FAILED
            call deallocate_table(table)
            return
        end if

        xc%initialized = .true.

        call deallocate_table(table)
    end subroutine xc_lsda_init

    !> Get exchange-correlation energy per particle at (n_up, n_dw)
    !!
    !! For U > 0 this is a direct evaluation of the tabulated functional
    !! (`get_exc_positive`). For U < 0 the Shiba partial particle-hole
    !! transformation is applied first:
    !!   e_xc(n_up, n_dw; U) = e_xc(1 - n_up, n_dw; |U|)
    !! matching `exc_value` in the C++ reference, where the recursion happens
    !! before any region/symmetry logic.
    !!
    !! @param[in]  xc   Initialized XC functional
    !! @param[in]  n_up Spin-up density (0 ≤ n_up ≤ 1)
    !! @param[in]  n_dw Spin-down density (0 ≤ n_dw ≤ 1)
    !! @param[out] exc  Exchange-correlation energy e_xc
    !! @param[out] ierr Error code (0 = success)
    subroutine get_exc(xc, n_up, n_dw, exc, ierr)
        type(xc_lsda_t), intent(in) :: xc
        real(dp), intent(in) :: n_up, n_dw
        real(dp), intent(out) :: exc
        integer, intent(out) :: ierr

        xc_evaluation_count = xc_evaluation_count + 1

        ! Check initialization
        if (.not. xc%initialized) then
            exc = 0.0_dp
            ierr = ERROR_INVALID_INPUT
            return
        end if

        ! Special case: U = 0 (free Fermi gas)
        if (abs(xc%U) < U_SMALL) then
            exc = 0.0_dp
            ierr = ERROR_SUCCESS
            return
        end if

        ! Check density bounds on the PHYSICAL densities, before any Shiba
        ! transformation (which would map an unphysical negative n_up back into
        ! the admissible window and hide the error).
        ! Note: Particle-hole symmetry allows n_up, n_down > 1.0
        ! Maximum total density is n_up + n_down ≤ 2 (Pauli exclusion)
        ! Use small tolerance to handle floating-point roundoff (e.g., 1.0+1.0 = 2.0000...009)
        if (n_up < -1.0e-12_dp .or. n_dw < -1.0e-12_dp .or. &
            (n_up + n_dw) > 2.0_dp + 1.0e-10_dp) then
            exc = 0.0_dp
            ierr = ERROR_OUT_OF_BOUNDS
            return
        end if

        if (xc%U < 0.0_dp) then
            ! Shiba: (n_up, n_dw; U < 0) -> (1 - n_up, n_dw; |U|), no sign flip
            call get_exc_positive(xc, 1.0_dp - n_up, n_dw, exc, ierr)
            return
        end if

        call get_exc_positive(xc, n_up, n_dw, exc, ierr)
    end subroutine get_exc

    !> Evaluate e_xc on the tabulated (repulsive) functional
    !!
    !! Thin wrapper around the recursive kernel `exc_recursive`; it only repeats
    !! the density bound check and converts the result into the (value, ierr)
    !! pair expected by `get_exc`. The Shiba transformation for U < 0 is NOT
    !! applied here: this routine always evaluates the |U| functional, so it
    !! must never be called by `get_exc` twice for the same point.
    !!
    !! @param[in]  xc   Initialized XC functional
    !! @param[in]  n_up Spin-up density (already Shiba-transformed if U < 0)
    !! @param[in]  n_dw Spin-down density
    !! @param[out] exc  Exchange-correlation energy e_xc
    !! @param[out] ierr Error code (0 = success)
    subroutine get_exc_positive(xc, n_up, n_dw, exc, ierr)
        type(xc_lsda_t), intent(in) :: xc
        real(dp), intent(in) :: n_up, n_dw
        real(dp), intent(out) :: exc
        integer, intent(out) :: ierr

        ! Check density bounds (see get_exc); still enforced here because the
        ! Shiba transformation can only produce admissible arguments from
        ! admissible input, so a failure means a genuine programming error.
        if (n_up < -1.0e-12_dp .or. n_dw < -1.0e-12_dp .or. &
            (n_up + n_dw) > 2.0_dp + 1.0e-10_dp) then
            exc = 0.0_dp
            ierr = ERROR_OUT_OF_BOUNDS
            return
        end if

        exc = exc_recursive(xc, n_up, n_dw)
        ierr = ERROR_SUCCESS
    end subroutine get_exc_positive

    !> Recursive kernel of e_xc for the tabulated (repulsive) functional
    !!
    !! Literal transcription of `exc_value` in the C++ reference
    !! (`original/spline2D.cc`): the empty-channel and corner shortcuts are
    !! evaluated at the TOP OF EVERY recursion level, so they see both the
    !! original densities and every intermediate pair produced by the region
    !! symmetries. That ordering is the whole point of using recursion here
    !! instead of one up-front call to `apply_symmetry_transform`: at n = 2 the
    !! point (1, 1) is in Region IV and only becomes an empty lattice AFTER the
    !! transform, so a shortcut tested on the original densities never fires and
    !! the spline is extrapolated into the empty corner instead of returning 0.
    !!
    !! In (n, m) variables the two C++ conditions read
    !!   |n - m| < 1e-14  <=>  n_dw = 0 (spin-down channel empty)
    !!   |n + m| < 1e-14  <=>  n_up = 0 (spin-up channel empty)
    !! and with one channel empty there is no double occupancy, hence e_xc = 0.
    !!
    !! The recursion is at most two levels deep: every region maps into Region I
    !! in a single step.
    !!
    !! @param[in] xc   Initialized XC functional
    !! @param[in] n_up Spin-up density
    !! @param[in] n_dw Spin-down density
    !! @return         e_xc(n_up, n_dw) for the |U| functional
    recursive function exc_recursive(xc, n_up_in, n_dw_in) result(exc)
        type(xc_lsda_t), intent(in) :: xc
        real(dp), intent(in) :: n_up_in, n_dw_in
        real(dp) :: exc

        integer :: region
        real(dp) :: n_up, n_dw, n_up_map, n_dw_map, n, m, n_up_c, n_dw_c

        ! Snap onto the region boundaries first, then use the snapped pair for
        ! the shortcuts, the region decision AND the evaluation (see
        ! `snap_to_boundaries`).
        call snap_to_boundaries(n_up_in, n_dw_in, n_up, n_dw, n, m)

        if (is_corner_shortcut(n, m)) then
            exc = 0.0_dp
            return
        end if

        ! Empty spin-down (|n - m| < tol) or empty spin-up (|n + m| < tol)
        if (abs(n - m) < EMPTY_CHANNEL_TOL .or. abs(n + m) < EMPTY_CHANNEL_TOL) then
            exc = 0.0_dp
            return
        end if

        region = determine_region(n, m)

        if (region /= 1) then
            ! Regions II, III and IV: recurse on the transformed pair, exactly as
            ! the C++ does. e_xc is invariant (no sign flip, no spin exchange).
            call apply_symmetry_transform(region, n_up, n_dw, n_up_map, n_dw_map)
            exc = exc_recursive(xc, n_up_map, n_dw_map)
            return
        end if

        ! Region I: clip to the tabulated range to absorb roundoff at the
        ! boundaries (the particle-hole transform can produce -1e-17)
        n_up_c = max(0.0_dp, min(1.0_dp, n_up))
        n_dw_c = max(0.0_dp, min(1.0_dp, n_dw))
        call convert_to_nm(n_up_c, n_dw_c, n, m)

        ! Fully polarized beyond the last tabulated density: e_xc = 0.
        ! Literal guard of the C++ `exc_value`; it also keeps the density
        ! spline below from degenerating to a zero-width window.
        if (m > xc%spl_exc%x(xc%spl_exc%n_x) - 1.0e-15_dp) then
            exc = 0.0_dp
            return
        end if

        ! Evaluate spline. The synthetic node at n = m is the fully polarized
        ! line, where e_xc vanishes; the prescribed derivative there is
        ! (∂e_xc/∂n_dn)/2 = (de_xc/dn) at fixed magnetization.
        exc = spline2d_eval(xc%spl_exc, n, m, &
                            node0_value = 0.0_dp, &
                            dfdx_first = 0.5_dp * dexc_dndown_b0(abs(xc%U), m))
    end function exc_recursive

    !> Get exchange-correlation potentials at (n_up, n_dw)
    !!
    !! Returns V_xc^up and V_xc^dn using symmetries.
    !!
    !! V_xc is genuinely discontinuous at n = n_up + n_dw = 1 (it is the
    !! derivative of the BALDA Mott gap): the branch used for n > 1 carries an
    !! overall minus sign, so the potential jumps by 2|v_base(1, m)| there
    !! (1.287 for U = 4, m = 0). The C++ reference has the same jump and the
    !! default configuration reproduces it exactly.
    !!
    !! If `xc%smoothing_width` (w) is positive, the jump is replaced, for
    !! |n - 1| < w, by a straight line between the values on the two sides of
    !! the discontinuity: V_xc is evaluated at n = 1 - w and n = 1 + w at fixed
    !! magnetization m and linearly interpolated in n. The result is continuous
    !! and coincides with the unsmoothed functional outside the window, at the
    !! price of departing from the C++ reference inside it.
    !!
    !! For U < 0 the Shiba partial particle-hole transformation is applied at
    !! the outermost level, before the region logic and before the smoothing:
    !!   V_xc^up(n_up, n_dw; U) = -V_xc^up(1 - n_up, n_dw; |U|)
    !!   V_xc^dn(n_up, n_dw; U) = +V_xc^dn(1 - n_up, n_dw; |U|)
    !! The asymmetry (up flips sign, down does not) is literal C++ behaviour
    !! (`Vxc_up_value` / `Vxc_dn_value` in `original/spline2D.cc`).
    !!
    !! @param[in]  xc      Initialized XC functional
    !! @param[in]  n_up    Spin-up density (0 ≤ n_up ≤ 1)
    !! @param[in]  n_dw    Spin-down density (0 ≤ n_dw ≤ 1)
    !! @param[out] v_xc_up XC potential for spin-up
    !! @param[out] v_xc_dw XC potential for spin-down
    !! @param[out] ierr    Error code (0 = success)
    subroutine get_vxc(xc, n_up, n_dw, v_xc_up, v_xc_dw, ierr)
        type(xc_lsda_t), intent(in) :: xc
        real(dp), intent(in) :: n_up, n_dw
        real(dp), intent(out) :: v_xc_up, v_xc_dw
        integer, intent(out) :: ierr

        real(dp) :: v_up_tmp, v_dw_tmp

        xc_evaluation_count = xc_evaluation_count + 1

        ! Check initialization
        if (.not. xc%initialized) then
            v_xc_up = 0.0_dp
            v_xc_dw = 0.0_dp
            ierr = ERROR_INVALID_INPUT
            return
        end if

        ! Special case: U = 0 (free Fermi gas)
        if (abs(xc%U) < U_SMALL) then
            v_xc_up = 0.0_dp
            v_xc_dw = 0.0_dp
            ierr = ERROR_SUCCESS
            return
        end if

        ! Check density bounds on the PHYSICAL densities, before any Shiba
        ! transformation (which would map an unphysical negative n_up back into
        ! the admissible window and hide the error).
        ! Note: Particle-hole symmetry allows n_up, n_down > 1.0
        ! Maximum total density is n_up + n_down ≤ 2 (Pauli exclusion)
        ! Use small tolerance to handle floating-point roundoff (e.g., 1.0+1.0 = 2.0000...009)
        if (n_up < -1.0e-12_dp .or. n_dw < -1.0e-12_dp .or. &
            (n_up + n_dw) > 2.0_dp + 1.0e-10_dp) then
            v_xc_up = 0.0_dp
            v_xc_dw = 0.0_dp
            ierr = ERROR_OUT_OF_BOUNDS
            return
        end if

        if (xc%U < 0.0_dp) then
            ! Shiba: (n_up, n_dw; U < 0) -> (1 - n_up, n_dw; |U|)
            call get_vxc_positive(xc, 1.0_dp - n_up, n_dw, v_up_tmp, v_dw_tmp, ierr)
            v_xc_up = -v_up_tmp
            v_xc_dw =  v_dw_tmp
            return
        end if

        call get_vxc_positive(xc, n_up, n_dw, v_xc_up, v_xc_dw, ierr)
    end subroutine get_vxc

    !> Evaluate V_xc on the tabulated (repulsive) functional
    !!
    !! Applies the region symmetries and the optional smoothing of the n = 1
    !! discontinuity of the |U| functional. The Shiba transformation for U < 0
    !! is NOT applied here; see `get_vxc`.
    !!
    !! @param[in]  xc      Initialized XC functional
    !! @param[in]  n_up    Spin-up density (already Shiba-transformed if U < 0)
    !! @param[in]  n_dw    Spin-down density
    !! @param[out] v_xc_up XC potential for spin-up
    !! @param[out] v_xc_dw XC potential for spin-down
    !! @param[out] ierr    Error code (0 = success)
    subroutine get_vxc_positive(xc, n_up, n_dw, v_xc_up, v_xc_dw, ierr)
        type(xc_lsda_t), intent(in) :: xc
        real(dp), intent(in) :: n_up, n_dw
        real(dp), intent(out) :: v_xc_up, v_xc_dw
        integer, intent(out) :: ierr

        real(dp) :: n, m, w, t
        real(dp) :: n_lo, n_hi, m_lo, m_hi
        real(dp) :: v_lo_up, v_lo_dw, v_hi_up, v_hi_dw

        ! Check density bounds (see get_vxc); still enforced here because the
        ! Shiba transformation can only produce admissible arguments from
        ! admissible input, so a failure means a genuine programming error.
        if (n_up < -1.0e-12_dp .or. n_dw < -1.0e-12_dp .or. &
            (n_up + n_dw) > 2.0_dp + 1.0e-10_dp) then
            v_xc_up = 0.0_dp
            v_xc_dw = 0.0_dp
            ierr = ERROR_OUT_OF_BOUNDS
            return
        end if

        ! NOTE: there is deliberately no "both densities zero" shortcut here.
        ! The empty lattice is already covered, per spin channel and exactly as
        ! in the C++ reference, by the shortcuts inside `vxc_up_recursive` /
        ! `vxc_dn_recursive` (and V_xc^dn(0, 0) = dexc_dndown_b0(|U|, 0) = 0).
        call convert_to_nm(n_up, n_dw, n, m)
        w = xc%smoothing_width

        if (w > 0.0_dp .and. abs(n - 1.0_dp) < w) then
            ! -----------------------------------------------------------------
            ! Smoothed branch: linear interpolation in n, at fixed m, between
            ! the two sides of the discontinuity.
            !
            ! The magnetization is clamped at each edge so that both mapped
            ! densities stay inside [0, 1]: |m| <= n at n = 1 - w and
            ! |m| <= 2 - n at n = 1 + w. This only bites for nearly fully
            ! polarized sites, where the window edge would otherwise leave the
            ! physical triangle.
            ! -----------------------------------------------------------------
            n_lo = 1.0_dp - w
            n_hi = 1.0_dp + w
            m_lo = max(-n_lo, min(n_lo, m))
            m_hi = max(-(2.0_dp - n_hi), min(2.0_dp - n_hi, m))

            call eval_vxc_branch(xc, 0.5_dp * (n_lo + m_lo), 0.5_dp * (n_lo - m_lo), v_lo_up, v_lo_dw)
            call eval_vxc_branch(xc, 0.5_dp * (n_hi + m_hi), 0.5_dp * (n_hi - m_hi), v_hi_up, v_hi_dw)

            t = (n - n_lo) / (2.0_dp * w)
            v_xc_up = (1.0_dp - t) * v_lo_up + t * v_hi_up
            v_xc_dw = (1.0_dp - t) * v_lo_dw + t * v_hi_dw
        else
            call eval_vxc_branch(xc, n_up, n_dw, v_xc_up, v_xc_dw)
        end if

        ierr = ERROR_SUCCESS
    end subroutine get_vxc_positive

    !> Evaluate V_xc at (n_up, n_dw) on the unsmoothed functional
    !!
    !! Thin wrapper that evaluates the two independent recursive kernels, one per
    !! spin channel, exactly as the C++ reference has two independent functions
    !! `Vxc_up_value` and `Vxc_dn_value`. They cannot share a single mapped point
    !! because their shortcuts differ (see `vxc_dn_recursive`).
    !!
    !! @param[in]  xc      Initialized XC functional
    !! @param[in]  n_up    Spin-up density
    !! @param[in]  n_dw    Spin-down density
    !! @param[out] v_xc_up XC potential for spin-up
    !! @param[out] v_xc_dw XC potential for spin-down
    subroutine eval_vxc_branch(xc, n_up, n_dw, v_xc_up, v_xc_dw)
        type(xc_lsda_t), intent(in) :: xc
        real(dp), intent(in) :: n_up, n_dw
        real(dp), intent(out) :: v_xc_up, v_xc_dw

        v_xc_up = vxc_up_recursive(xc, n_up, n_dw)
        v_xc_dw = vxc_dn_recursive(xc, n_up, n_dw)
    end subroutine eval_vxc_branch

    !> Fully polarized corner shortcut of the C++ reference
    !!
    !! Reproduces, literally,
    !!   `if((1.0-1.0e-14 <= dens && dens <= 1.0) && (1.0-1.0e-14 <= mag && mag <= 1.0)) return 0.0;`
    !! which appears at the top of all three C++ evaluators (`exc_value`,
    !! `Vxc_up_value`, `Vxc_dn_value`), BEFORE the empty-channel shortcut. It
    !! catches the fully polarized half-filled point n = m = 1 (n_up = 1,
    !! n_dw = 0) in a band that the empty-channel test, which measures
    !! |n - m| = 2 n_dw, would just miss.
    !!
    !! @param[in] n Total density n = n_up + n_dw
    !! @param[in] m Magnetization m = n_up - n_dw
    !! @return      .true. if the shortcut applies (value is 0)
    function is_corner_shortcut(n, m) result(hit)
        real(dp), intent(in) :: n, m
        logical :: hit

        hit = (1.0_dp - CORNER_TOL <= n .and. n <= 1.0_dp) .and. &
              (1.0_dp - CORNER_TOL <= m .and. m <= 1.0_dp)
    end function is_corner_shortcut

    !> Recursive kernel of V_xc^up for the tabulated (repulsive) functional
    !!
    !! Literal transcription of `Vxc_up_value` in the C++ reference. As in
    !! `exc_recursive`, the shortcuts are re-tested at every recursion level, so
    !! they also see the pairs produced by the region symmetries (that is how
    !! (1, 1) in Region IV reaches the empty lattice (0, 0) and returns 0).
    !!
    !! Note the asymmetry with respect to `exc_recursive`: only the empty
    !! spin-DOWN test (|n - m| < tol) short-circuits here, the C++ has no
    !! `|n + m|` clause in this function.
    !!
    !! @param[in] xc   Initialized XC functional
    !! @param[in] n_up Spin-up density
    !! @param[in] n_dw Spin-down density
    !! @return         V_xc^up(n_up, n_dw) for the |U| functional
    recursive function vxc_up_recursive(xc, n_up_in, n_dw_in) result(v_up)
        type(xc_lsda_t), intent(in) :: xc
        real(dp), intent(in) :: n_up_in, n_dw_in
        real(dp) :: v_up

        integer :: region
        real(dp) :: n_up, n_dw, n_up_map, n_dw_map, n, m, n_up_c, n_dw_c

        ! Snap onto the region boundaries first (see `snap_to_boundaries`)
        call snap_to_boundaries(n_up_in, n_dw_in, n_up, n_dw, n, m)

        if (is_corner_shortcut(n, m)) then
            v_up = 0.0_dp
            return
        end if

        ! Empty spin-down channel: no double occupancy, V_xc^up = 0
        if (abs(n - m) < EMPTY_CHANNEL_TOL) then
            v_up = 0.0_dp
            return
        end if

        region = determine_region(n, m)

        if (region /= 1) then
            call apply_symmetry_transform(region, n_up, n_dw, n_up_map, n_dw_map)
            select case (region)
            case (2)
                ! Region II: spin exchange turns V_xc^up into the SIBLING function
                v_up = vxc_dn_recursive(xc, n_up_map, n_dw_map)
            case (3)
                ! Region III: particle-hole, same spin, overall sign flip
                v_up = -vxc_up_recursive(xc, n_up_map, n_dw_map)
            case default
                ! Region IV: particle-hole + spin exchange, sign flip
                v_up = -vxc_dn_recursive(xc, n_up_map, n_dw_map)
            end select
            return
        end if

        ! Region I: clip to the tabulated range to absorb boundary roundoff
        n_up_c = max(0.0_dp, min(1.0_dp, n_up))
        n_dw_c = max(0.0_dp, min(1.0_dp, n_dw))
        call convert_to_nm(n_up_c, n_dw_c, n, m)

        ! The three asymmetries with `vxc_dn_recursive` are deliberate parity
        ! with the C++ reference (`original/spline2D.cc`), NOT oversights:
        !  - the "fully polarized beyond the table" guard below exists in
        !    `Vxc_up_value` but is ABSENT from `Vxc_dn_value`;
        !  - the linear fallback for a single-row window exists in
        !    `Vxc_up_value` but is ABSENT from `Vxc_dn_value`;
        !  - the synthetic node at n = m vanishes for V_xc^up (e_xc is
        !    identically zero along n_dn = 0) but equals ∂e_xc/∂n_dn for
        !    V_xc^dn, which is the analytic `dexc_dndown_b0`.
        if (m > xc%spl_vxc_up%x(xc%spl_vxc_up%n_x) - 1.0e-15_dp) then
            v_up = 0.0_dp
            return
        end if

        v_up = spline2d_eval(xc%spl_vxc_up, n, m, node0_value = 0.0_dp)
    end function vxc_up_recursive

    !> Recursive kernel of V_xc^dn for the tabulated (repulsive) functional
    !!
    !! Literal transcription of `Vxc_dn_value` in the C++ reference.
    !!
    !! DELIBERATE DIVERGENCE FROM THE TASK TEXT (T18 of NEXT_STEPS_REPORT.md):
    !! the report asks for "zero for exc, v_xc_up and v_xc_dw" when a channel is
    !! empty. That generalization is wrong for the spin-down potential. With the
    !! spin-down channel empty the C++ returns
    !!   `return (dexc_dndownB0(u, mag));`   (`original/spline2D.cc:677`)
    !! and NOT zero: adding one spin-down electron to a polarized band does cost
    !! correlation energy, so ∂e_xc/∂n_dn is finite there even though e_xc itself
    !! and V_xc^up vanish. Returning zero would truncate the only XC force that
    !! repopulates an empty channel. The value happens to vanish at m = 0, which
    !! is why V_xc^dn(1, 1) = 0 all the same (Region IV -> -V_xc^up(0, 0) = 0).
    !!
    !! There is also no `|n + m|` (empty spin-up) clause in the C++ here, and no
    !! "beyond the last tabulated density" guard; both absences are reproduced.
    !!
    !! @param[in] xc   Initialized XC functional
    !! @param[in] n_up Spin-up density
    !! @param[in] n_dw Spin-down density
    !! @return         V_xc^dn(n_up, n_dw) for the |U| functional
    recursive function vxc_dn_recursive(xc, n_up_in, n_dw_in) result(v_dw)
        type(xc_lsda_t), intent(in) :: xc
        real(dp), intent(in) :: n_up_in, n_dw_in
        real(dp) :: v_dw

        integer :: region
        real(dp) :: n_up, n_dw, n_up_map, n_dw_map, n, m, n_up_c, n_dw_c

        ! Snap onto the region boundaries first (see `snap_to_boundaries`)
        call snap_to_boundaries(n_up_in, n_dw_in, n_up, n_dw, n, m)

        if (is_corner_shortcut(n, m)) then
            v_dw = 0.0_dp
            return
        end if

        ! Empty spin-down channel: at n_dn = 0 one has n = m, and V_xc^dn is the
        ! analytic ∂e_xc/∂n_dn of the fully polarized band (see above).
        if (abs(n - m) < EMPTY_CHANNEL_TOL) then
            v_dw = dexc_dndown_b0(abs(xc%U), m)
            return
        end if

        region = determine_region(n, m)

        if (region /= 1) then
            call apply_symmetry_transform(region, n_up, n_dw, n_up_map, n_dw_map)
            select case (region)
            case (2)
                ! Region II: spin exchange turns V_xc^dn into the SIBLING function
                v_dw = vxc_up_recursive(xc, n_up_map, n_dw_map)
            case (3)
                v_dw = -vxc_dn_recursive(xc, n_up_map, n_dw_map)
            case default
                v_dw = -vxc_up_recursive(xc, n_up_map, n_dw_map)
            end select
            return
        end if

        ! Region I: clip to the tabulated range to absorb boundary roundoff
        n_up_c = max(0.0_dp, min(1.0_dp, n_up))
        n_dw_c = max(0.0_dp, min(1.0_dp, n_dw))
        call convert_to_nm(n_up_c, n_dw_c, n, m)

        v_dw = spline2d_eval(xc%spl_vxc_down, n, m, &
                             node0_value = dexc_dndown_b0(abs(xc%U), m), &
                             allow_linear_branch = .false.)
    end function vxc_dn_recursive

    !> Destroy XC functional and free memory
    !!
    !! @param[inout] xc XC functional object
    subroutine xc_lsda_destroy(xc)
        type(xc_lsda_t), intent(inout) :: xc

        call spline2d_destroy(xc%spl_exc)
        call spline2d_destroy(xc%spl_vxc_up)
        call spline2d_destroy(xc%spl_vxc_down)

        xc%initialized = .false.
        xc%U = 0.0_dp
        xc%smoothing_width = 0.0_dp
    end subroutine xc_lsda_destroy

    !> Auxiliary Bethe Ansatz integral of the fully polarized limit
    !!
    !! Computes, in closed form,
    !!
    !!   F(γ, Q) = 2 ∫₀^Q γ cos²k / (γ² + sin²k) dk
    !!           = 2 [ √(1+γ²) · A(Q) - γ Q ],
    !!   A(Q) = arctan( √(1+γ²) tan Q / γ )   (continued across Q = π/2)
    !!
    !! This is the function `integral_1(gama, Q)` of the C++ reference
    !! (`original/spline2D.cc`), used only by `dexc_dndown_b0`. The closed form
    !! was verified against a high-order quadrature of the integral above and
    !! against the fully polarized column of the tabulated V_xc^dn.
    !!
    !! @param[in] gama Interaction parameter γ = |U|/4, strictly positive
    !! @param[in] q    Upper limit of integration, Q = π n ∈ [0, π]
    !! @return         Value of F(γ, Q)
    function integral_1(gama, q) result(res)
        real(dp), intent(in) :: gama, q
        real(dp) :: res

        real(dp) :: s, a

        s = sqrt(1.0_dp + gama * gama)

        if (abs(q - 0.5_dp * PI) < 1.0e-14_dp) then
            ! tan(Q) diverges; A(π/2) = π/2 by continuity
            a = 0.5_dp * PI
        else if (q < 0.5_dp * PI) then
            a = atan(s * tan(q) / gama)
        else
            ! Second branch of the arctangent: A must stay monotonic in Q
            a = atan(s * tan(q) / gama) + PI
        end if

        res = 2.0_dp * (s * a - gama * q)
    end function integral_1

    !> Derivative ∂e_xc/∂n_dn at the fully polarized point (n_up = n, n_dn = 0)
    !!
    !! Analytical Bethe Ansatz result for the cost of adding one spin-down
    !! electron to a fully polarized band of density n, with the Hartree term
    !! U·n already subtracted (the tables are stored the same way):
    !!
    !!   ∂e_xc/∂n_dn = 2 cos Q [2 arctan(sin Q / γ)/π - 1]
    !!                 + 2 [1 - F(γ, Q)/π] - U n,      Q = π n, γ = U/4
    !!
    !! This is `dexc_dndownB0(u, dens)` of the C++ reference. It supplies the
    !! value of the synthetic node of the density spline of V_xc^dn and, halved,
    !! the prescribed first derivative of the density spline of e_xc (see
    !! `spline2d_eval`).
    !!
    !! @param[in] u    Hubbard interaction, must be positive (|U| of the table)
    !! @param[in] dens Density of the polarized channel, n ∈ [0, 1]
    !! @return         ∂e_xc/∂n_dn at (n_up = dens, n_dn = 0); 0 if u ≈ 0
    function dexc_dndown_b0(u, dens) result(res)
        real(dp), intent(in) :: u, dens
        real(dp) :: res

        real(dp) :: gama, q, f

        if (abs(u) < U_SMALL) then
            ! e_xc ≡ 0 for the non-interacting system
            res = 0.0_dp
            return
        end if

        gama = u / 4.0_dp
        q = PI * dens
        f = integral_1(gama, q)

        res = 2.0_dp * cos(q) * (2.0_dp * atan(sin(q) / gama) / PI - 1.0_dp) &
              + 2.0_dp * (1.0_dp - f / PI) - u * dens
    end function dexc_dndown_b0

    !> Determine symmetry region for a point given in (n, m) coordinates
    !!
    !! Region I:   m ≥ 0 and n ≤ 1
    !! Region II:  m < 0 and n ≤ 1
    !! Region III: m < 0 and n > 1
    !! Region IV:  m ≥ 0 and n > 1
    !!
    !! The comparisons are STRICT, exactly as in the C++ reference
    !! (`original/spline2D.cc`, tests `mag < 0.0 && dens <= 1.0` etc.). Roundoff
    !! at the boundaries must be absorbed BEFORE this function is called, by
    !! `snap_to_boundaries`, and the snapped (n, m) must then also be the pair
    !! that is fed to the spline: deciding the branch with a tolerance while
    !! evaluating at the untolerated point is what turned a particle-hole
    !! transformation into a spin exchange (and vice versa) in a 1e-12 band
    !! around n = 1 and m = 0.
    !!
    !! @param[in] n Total density n = n_up + n_dw (already snapped)
    !! @param[in] m Magnetization m = n_up - n_dw (already snapped)
    !! @return      Region number (1-4)
    function determine_region(n, m) result(region)
        real(dp), intent(in) :: n, m
        integer :: region

        if (m >= 0.0_dp .and. n <= 1.0_dp) then
            region = 1
        else if (m < 0.0_dp .and. n <= 1.0_dp) then
            region = 2
        else if (m < 0.0_dp) then
            region = 3  ! m < 0 and n > 1
        else
            region = 4  ! m ≥ 0 and n > 1
        end if
    end function determine_region

    !> Snap a density pair onto the region boundaries n = 1 and m = 0
    !!
    !! Returns the (n, m) pair to be used BOTH for the region decision and for
    !! the spline evaluation, together with the matching (n_up, n_dw).
    !! Whenever |n - 1| or |m| is below `REGION_SNAP_TOL`, or |m| exceeds n by
    !! less than `REGION_SNAP_TOL`, the coordinate is fixed at the boundary
    !! value and the spin densities are recomputed from the snapped pair, so
    !! that region and evaluation can never disagree.
    !!
    !! Why the snap is MANDATORY (not cosmetic): the four region branches are
    !! different functions, not neighbouring values of one function. At
    !! n = 1 + eps with m < 0, Region II returns V_xc of the OTHER spin channel
    !! with no sign change while Region III returns the SAME channel with the
    !! sign flipped. Choosing the region with a tolerance (so that n = 1 + eps
    !! counts as n ≤ 1) while evaluating the spline at the raw n = 1 + eps
    !! therefore produced a full-size (order 1.66 for |U| = 4) potential in the
    !! wrong spin channel and with the wrong sign for every eps in
    !! (0, REGION_SNAP_TOL]. Snapping n to exactly 1 makes both the region test
    !! and the evaluation see the same point, so the answer is the n = 1 limit
    !! of the branch actually used. The same argument applies to m = 0, where
    !! the tolerance would otherwise break the identity
    !! V_xc^up(n_up, n_dw) = V_xc^dn(n_dw, n_up).
    !!
    !! The spin densities are only rebuilt when the point sits inside one of the
    !! snap bands (a point already exactly on the boundary is rebuilt too, which
    !! is a no-op up to one ulp); points away from the boundaries are passed
    !! through bit for bit.
    !!
    !! @param[in]  n_up_in Raw spin-up density
    !! @param[in]  n_dw_in Raw spin-down density
    !! @param[out] n_up    Spin-up density consistent with the snapped (n, m)
    !! @param[out] n_dw    Spin-down density consistent with the snapped (n, m)
    !! @param[out] n       Snapped total density
    !! @param[out] m       Snapped magnetization
    subroutine snap_to_boundaries(n_up_in, n_dw_in, n_up, n_dw, n, m)
        real(dp), intent(in) :: n_up_in, n_dw_in
        real(dp), intent(out) :: n_up, n_dw, n, m

        logical :: snapped

        call convert_to_nm(n_up_in, n_dw_in, n, m)

        snapped = .false.

        if (abs(n - 1.0_dp) < REGION_SNAP_TOL) then
            n = 1.0_dp
            snapped = .true.
        end if

        if (abs(m) < REGION_SNAP_TOL) then
            m = 0.0_dp
            snapped = .true.
        end if

        ! Third boundary: the physical triangle |m| <= n (one spin density cannot
        ! be negative). Roundoff can push a site to n_up = 1 + 4e-16 at the full
        ! band, which after the Shiba transformation becomes n_dw = -2e-16 and
        ! hence |m| = n + 4e-16. That tiny excursion outside the triangle skips
        ! the fully polarized corner shortcut (whose test is `m <= 1`) and lands
        ! on the empty-channel shortcut instead, which returns the full-size
        ! dexc_dndown_b0 (1.66 for |U| = 4) where 0 is due. Snapping |m| back to
        ! n keeps the evaluation inside the tabulated domain 0 <= m <= n.
        if (abs(m) > n .and. abs(m) - n < REGION_SNAP_TOL) then
            m = sign(n, m)
            snapped = .true.
        end if

        if (snapped) then
            ! (n, m) are kept as the exact snapped values; the spin densities are
            ! rebuilt from them. Recomputing n = n_up + n_dw afterwards would
            ! re-introduce the very 1-ulp drift that was just removed, so it is
            ! deliberately not done.
            n_up = 0.5_dp * (n + m)
            n_dw = 0.5_dp * (n - m)
        else
            n_up = n_up_in
            n_dw = n_dw_in
        end if
    end subroutine snap_to_boundaries

    !> Convert (n_up, n_dw) to (n, m) coordinates
    !!
    !! n = n_up + n_dw (total density)
    !! m = n_up - n_dw (magnetization)
    !!
    !! @param[in]  n_up Spin-up density
    !! @param[in]  n_dw Spin-down density
    !! @param[out] n    Total density
    !! @param[out] m    Magnetization
    subroutine convert_to_nm(n_up, n_dw, n, m)
        real(dp), intent(in) :: n_up, n_dw
        real(dp), intent(out) :: n, m

        n = n_up + n_dw
        m = n_up - n_dw
    end subroutine convert_to_nm

    !> Apply symmetry transformation to map any region to Region I
    !!
    !! Transformations:
    !! Region I:   (n_up, n_dw) → (n_up, n_dw)
    !! Region II:  (n_up, n_dw) → (n_dw, n_up)     [spin exchange]
    !! Region III: (n_up, n_dw) → (1-n_up, 1-n_dw) [particle-hole]
    !! Region IV:  (n_up, n_dw) → (1-n_dw, 1-n_up) [combined]
    !!
    !! @param[in]  region   Symmetry region (1-4)
    !! @param[in]  n_up     Original spin-up density
    !! @param[in]  n_dw     Original spin-down density
    !! @param[out] n_up_map Mapped spin-up density
    !! @param[out] n_dw_map Mapped spin-down density
    subroutine apply_symmetry_transform(region, n_up, n_dw, n_up_map, n_dw_map)
        integer, intent(in) :: region
        real(dp), intent(in) :: n_up, n_dw
        real(dp), intent(out) :: n_up_map, n_dw_map

        select case (region)
        case (1)
            ! Region I: Direct (no transformation)
            n_up_map = n_up
            n_dw_map = n_dw

        case (2)
            ! Region II: Spin exchange
            n_up_map = n_dw
            n_dw_map = n_up

        case (3)
            ! Region III: Particle-hole symmetry
            n_up_map = 1.0_dp - n_up
            n_dw_map = 1.0_dp - n_dw

        case (4)
            ! Region IV: Combined (particle-hole + spin exchange)
            n_up_map = 1.0_dp - n_dw
            n_dw_map = 1.0_dp - n_up

        case default
            ! Should never happen (region is always 1-4 from determine_region)
            ! Return identity transformation as fallback
            n_up_map = n_up
            n_dw_map = n_dw
        end select
    end subroutine apply_symmetry_transform
end module xc_lsda
