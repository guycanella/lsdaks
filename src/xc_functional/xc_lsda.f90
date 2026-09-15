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
module xc_lsda
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use lsda_constants, only: dp, U_SMALL
    use spline2d, only: spline2d_t, spline2d_init, spline2d_eval, spline2d_destroy
    use table_io, only: xc_table_t, read_fortran_table, deallocate_table
    use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT, ERROR_OUT_OF_BOUNDS, ERROR_FILE_READ, ERROR_SPLINE_INITIALIZATION_FAILED
    implicit none
    private

    !> LSDA XC functional type containing splines for exc, vxc_up, vxc_down
    type, public :: xc_lsda_t
        real(dp) :: U = 0.0_dp                  !< Hubbard U parameter
        real(dp) :: smoothing_width = 0.0_dp    !< Half-width w of the linear smoothing of V_xc around n = 1 (0 = off)
        type(spline2d_t) :: spl_exc             !< Spline for e_xc(n, m)
        type(spline2d_t) :: spl_vxc_up          !< Spline for V_xc^up(n, m)
        type(spline2d_t) :: spl_vxc_down        !< Spline for V_xc^dn(n, m)
        logical :: initialized = .false.        !< Initialization flag
    end type xc_lsda_t

    !> Largest admissible smoothing half-width (the window must stay inside 0 < n < 2)
    real(dp), parameter, public :: XC_SMOOTHING_WIDTH_MAX = 1.0_dp

    public :: xc_lsda_init
    public :: get_exc
    public :: get_vxc
    public :: xc_lsda_destroy

    private :: determine_region
    private :: convert_to_nm
    private :: apply_symmetry_transform
    private :: eval_vxc_branch

contains

    !> Initialize XC functional from table file
    !!
    !! Loads table and constructs 2D splines for exc, vxc_up, vxc_down.
    !!
    !! The optional `smoothing_width` activates the linear smoothing of the
    !! V_xc discontinuity at n = 1 (see `get_vxc`). The default, 0, disables it
    !! and reproduces the C++ reference exactly.
    !!
    !! @param[out] xc         XC functional object
    !! @param[in]  table_file Path to table file (Fortran binary format)
    !! @param[out] ierr Error code (0 = success)
    !! @param[in]  smoothing_width Optional half-width w of the V_xc smoothing
    !!                             window around n = 1; must be a non-NaN value
    !!                             with 0 <= w < XC_SMOOTHING_WIDTH_MAX
    !!                             (default 0 = off). NaN is rejected with
    !!                             ERROR_INVALID_INPUT because it would otherwise
    !!                             disable the smoothing silently in `get_vxc`.
    subroutine xc_lsda_init(xc, table_file, ierr, smoothing_width)
        type(xc_lsda_t), intent(out) :: xc
        character(len=*), intent(in) :: table_file
        integer, intent(out) :: ierr
        real(dp), intent(in), optional :: smoothing_width

        type(xc_table_t) :: table
        integer :: io_stat
        integer, allocatable :: n_y_pts(:)
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

        call read_fortran_table(table_file, table, io_stat)
        if (io_stat /= 0) then
            print *, "ERROR: Failed to read table file: ", trim(table_file)
            ierr = ERROR_FILE_READ
            return
        end if

        xc%U = table%U

        ! Create n_y_pts array (same for all n_grid points in current implementation)
        allocate(n_y_pts(table%n_points_n))
        do i = 1, table%n_points_n
            n_y_pts(i) = table%n_points_m
        end do

        call spline2d_init(xc%spl_exc, table%n_grid, table%m_grid, &
                           table%exc, n_y_pts)

        call spline2d_init(xc%spl_vxc_up, table%n_grid, table%m_grid, &
                           table%vxc_up, n_y_pts)

        call spline2d_init(xc%spl_vxc_down, table%n_grid, table%m_grid, &
                           table%vxc_down, n_y_pts)

        deallocate(n_y_pts)

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
    !! Uses symmetries to map any (n_up, n_dw) to Region I and interpolates.
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

        integer :: region
        real(dp) :: n_up_map, n_dw_map, n, m

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

        ! Check density bounds
        ! Note: Particle-hole symmetry allows n_up, n_down > 1.0
        ! Maximum total density is n_up + n_down ≤ 2 (Pauli exclusion)
        ! Use small tolerance to handle floating-point roundoff (e.g., 1.0+1.0 = 2.0000...009)
        if (n_up < -1.0e-12_dp .or. n_dw < -1.0e-12_dp .or. &
            (n_up + n_dw) > 2.0_dp + 1.0e-10_dp) then
            exc = 0.0_dp
            ierr = ERROR_OUT_OF_BOUNDS
            return
        end if

        ! Special case: both densities zero
        if (n_up < 1.0e-12_dp .and. n_dw < 1.0e-12_dp) then
            exc = 0.0_dp
            ierr = ERROR_SUCCESS
            return
        end if

        ! Determine region and apply symmetry
        region = determine_region(n_up, n_dw)
        call apply_symmetry_transform(region, n_up, n_dw, n_up_map, n_dw_map)

        ! Clip to valid range to handle numerical errors near boundaries
        ! This is especially important for n ≈ 2.0 where particle-hole transform
        ! can produce small negative values due to roundoff
        n_up_map = max(0.0_dp, min(1.0_dp, n_up_map))
        n_dw_map = max(0.0_dp, min(1.0_dp, n_dw_map))

        call convert_to_nm(n_up_map, n_dw_map, n, m)

        ! Evaluate spline
        exc = spline2d_eval(xc%spl_exc, n, m)
        ierr = ERROR_SUCCESS
    end subroutine get_exc

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

        real(dp) :: n, m, w, t
        real(dp) :: n_lo, n_hi, m_lo, m_hi
        real(dp) :: v_lo_up, v_lo_dw, v_hi_up, v_hi_dw

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

        ! Check density bounds
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

        ! Special case: both densities zero
        if (n_up < 1.0e-12_dp .and. n_dw < 1.0e-12_dp) then
            v_xc_up = 0.0_dp
            v_xc_dw = 0.0_dp
            ierr = ERROR_SUCCESS
            return
        end if

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
    end subroutine get_vxc

    !> Evaluate V_xc at (n_up, n_dw) on the unsmoothed functional
    !!
    !! Maps the point to Region I with the physical symmetries, evaluates the
    !! tabulated splines and applies the region-specific sign/spin exchange.
    !! This is the exact C++ behaviour and carries the discontinuity at n = 1.
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

        integer :: region
        real(dp) :: n_up_map, n_dw_map, n, m
        real(dp) :: v_up_base, v_dw_base

        ! Determine region and apply symmetry
        region = determine_region(n_up, n_dw)
        call apply_symmetry_transform(region, n_up, n_dw, n_up_map, n_dw_map)

        ! Clip to valid range to handle numerical errors near boundaries
        n_up_map = max(0.0_dp, min(1.0_dp, n_up_map))
        n_dw_map = max(0.0_dp, min(1.0_dp, n_dw_map))

        call convert_to_nm(n_up_map, n_dw_map, n, m)

        ! Evaluate splines
        v_up_base = spline2d_eval(xc%spl_vxc_up, n, m)
        v_dw_base = spline2d_eval(xc%spl_vxc_down, n, m)

        ! Apply region-specific transformations
        select case (region)
        case (1)
            ! Region I (m ≥ 0, n ≤ 1): Direct
            v_xc_up = v_up_base
            v_xc_dw = v_dw_base

        case (2)
            ! Region II (m < 0, n ≤ 1): Spin exchange
            v_xc_up = v_dw_base
            v_xc_dw = v_up_base

        case (3)
            ! Region III (m < 0, n > 1): Particle-hole symmetry
            v_xc_up = -v_up_base
            v_xc_dw = -v_dw_base

        case (4)
            ! Region IV (m ≥ 0, n > 1): Combined
            v_xc_up = -v_dw_base
            v_xc_dw = -v_up_base

        case default
            v_xc_up = 0.0_dp
            v_xc_dw = 0.0_dp
        end select
    end subroutine eval_vxc_branch

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

    !> Determine symmetry region for (n_up, n_dw)
    !!
    !! Region I:   m ≥ 0 and n ≤ 1
    !! Region II:  m < 0 and n ≤ 1
    !! Region III: m < 0 and n > 1
    !! Region IV:  m ≥ 0 and n > 1
    !!
    !! @param[in] n_up Spin-up density
    !! @param[in] n_dw Spin-down density
    !! @return         Region number (1-4)
    function determine_region(n_up, n_dw) result(region)
        real(dp), intent(in) :: n_up, n_dw
        integer :: region

        real(dp) :: n, m
        real(dp), parameter :: TOL = 1.0e-12_dp  ! Tolerance for boundary cases

        n = n_up + n_dw
        m = n_up - n_dw

        ! Use tolerances to handle floating-point errors at n = 1.0 boundary
        ! This is critical for half-filling cases where n = n_up + n_down = 1.0
        if (m >= -TOL .and. n <= 1.0_dp + TOL) then
            region = 1
        else if (m < -TOL .and. n <= 1.0_dp + TOL) then
            region = 2
        else if (m < -TOL .and. n > 1.0_dp + TOL) then
            region = 3
        else
            region = 4  ! m ≥ 0 and n > 1
        end if
    end function determine_region

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
