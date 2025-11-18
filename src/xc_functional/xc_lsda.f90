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
    use lsda_constants, only: dp, U_SMALL
    use spline2d, only: spline2d_t, spline2d_init, spline2d_eval, spline2d_destroy
    use table_io, only: xc_table_t, read_fortran_table, deallocate_table
    use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT, ERROR_OUT_OF_BOUNDS, ERROR_FILE_READ, ERROR_SPLINE_INITIALIZATION_FAILED
    implicit none
    private

    !> LSDA XC functional type containing splines for exc, vxc_up, vxc_dw
    type, public :: xc_lsda_t
        real(dp) :: U = 0.0_dp                  !< Hubbard U parameter
        type(spline2d_t) :: spl_exc             !< Spline for e_xc(n, m)
        type(spline2d_t) :: spl_vxc_up          !< Spline for V_xc^up(n, m)
        type(spline2d_t) :: spl_vxc_dw          !< Spline for V_xc^dn(n, m)
        logical :: initialized = .false.        !< Initialization flag
    end type xc_lsda_t

    public :: xc_lsda_init
    public :: get_exc
    public :: get_vxc
    public :: xc_lsda_destroy

    private :: determine_region
    private :: convert_to_nm
    private :: apply_symmetry_transform

contains

    !> Initialize XC functional from table file
    !!
    !! Loads table and constructs 2D splines for exc, vxc_up, vxc_dw.
    !!
    !! @param[out] xc         XC functional object
    !! @param[in]  table_file Path to table file (Fortran binary format)
    !! @param[out] ierr Error code (0 = success)
    subroutine xc_lsda_init(xc, table_file, ierr)
        type(xc_lsda_t), intent(out) :: xc
        character(len=*), intent(in) :: table_file
        integer, intent(out) :: ierr

        type(xc_table_t) :: table
        integer :: io_stat
        integer, allocatable :: n_y_pts(:)
        integer :: i

        ierr = ERROR_SUCCESS

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

        call spline2d_init(xc%spl_vxc_dw, table%n_grid, table%m_grid, &
                           table%vxc_dw, n_y_pts)

        deallocate(n_y_pts)

        if (.not. xc%spl_exc%initialized .or. &
            .not. xc%spl_vxc_up%initialized .or. &
            .not. xc%spl_vxc_dw%initialized) then
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
        if (n_up < 0.0_dp .or. n_up > 1.0_dp .or. &
            n_dw < 0.0_dp .or. n_dw > 1.0_dp) then
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
        call convert_to_nm(n_up_map, n_dw_map, n, m)

        ! Evaluate spline
        exc = spline2d_eval(xc%spl_exc, n, m)
        ierr = ERROR_SUCCESS
    end subroutine get_exc

    !> Get exchange-correlation potentials at (n_up, n_dw)
    !!
    !! Returns V_xc^up and V_xc^dn using symmetries.
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

        integer :: region
        real(dp) :: n_up_map, n_dw_map, n, m
        real(dp) :: v_up_base, v_dw_base

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
        if (n_up < 0.0_dp .or. n_up > 1.0_dp .or. &
            n_dw < 0.0_dp .or. n_dw > 1.0_dp) then
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

        ! Determine region and apply symmetry
        region = determine_region(n_up, n_dw)
        call apply_symmetry_transform(region, n_up, n_dw, n_up_map, n_dw_map)
        call convert_to_nm(n_up_map, n_dw_map, n, m)

        ! Evaluate splines
        v_up_base = spline2d_eval(xc%spl_vxc_up, n, m)
        v_dw_base = spline2d_eval(xc%spl_vxc_dw, n, m)

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
        end select

        ierr = ERROR_SUCCESS
    end subroutine get_vxc

    !> Destroy XC functional and free memory
    !!
    !! @param[inout] xc XC functional object
    subroutine xc_lsda_destroy(xc)
        type(xc_lsda_t), intent(inout) :: xc

        call spline2d_destroy(xc%spl_exc)
        call spline2d_destroy(xc%spl_vxc_up)
        call spline2d_destroy(xc%spl_vxc_dw)

        xc%initialized = .false.
        xc%U = 0.0_dp
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

        n = n_up + n_dw
        m = n_up - n_dw

        if (m >= 0.0_dp .and. n <= 1.0_dp) then
            region = 1
        else if (m < 0.0_dp .and. n <= 1.0_dp) then
            region = 2
        else if (m < 0.0_dp .and. n > 1.0_dp) then
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
