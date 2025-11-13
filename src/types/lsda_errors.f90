!> Module for centralized error handling
!!
!! Provides consistent error codes, messages, and handling across the project.
!! All modules should use these error codes for uniformity.
!!
!! Error code ranges:
!! - 0: Success
!! - 1-99: Input validation errors
!! - 100-199: Numerical/computational errors
!! - 200-299: I/O errors
!! - 300-399: Memory allocation errors
module lsda_errors
    use lsda_constants, only: dp
    implicit none
    private

    ! Public error handling routines
    public :: error_handler
    public :: check_bounds
    public :: check_positive
    public :: check_range
    public :: get_error_message

    ! Public error codes
    public :: ERROR_SUCCESS
    public :: ERROR_INVALID_INPUT
    public :: ERROR_OUT_OF_BOUNDS
    public :: ERROR_SIZE_MISMATCH
    public :: ERROR_INVALID_CONCENTRATION
    public :: ERROR_NEGATIVE_VALUE
    public :: ERROR_INVALID_RANGE
    public :: ERROR_CONVERGENCE_FAILED
    public :: ERROR_SINGULAR_MATRIX
    public :: ERROR_LAPACK_FAILED
    public :: ERROR_FILE_NOT_FOUND
    public :: ERROR_FILE_READ
    public :: ERROR_FILE_WRITE
    public :: ERROR_ALLOCATION_FAILED
    public :: ERROR_DEALLOCATION_FAILED

    ! Error codes - General
    integer, parameter :: ERROR_SUCCESS = 0
    integer, parameter :: ERROR_INVALID_INPUT = 1
    integer, parameter :: ERROR_OUT_OF_BOUNDS = 2
    integer, parameter :: ERROR_SIZE_MISMATCH = 3
    integer, parameter :: ERROR_INVALID_CONCENTRATION = 4
    integer, parameter :: ERROR_NEGATIVE_VALUE = 5
    integer, parameter :: ERROR_INVALID_RANGE = 6
    integer, parameter :: ERROR_NOT_A_NUMBER = 7

    ! Error codes - Numerical (100-199)
    integer, parameter :: ERROR_CONVERGENCE_FAILED = 100
    integer, parameter :: ERROR_SINGULAR_MATRIX = 101
    integer, parameter :: ERROR_LAPACK_FAILED = 102

    ! Error codes - I/O (200-299)
    integer, parameter :: ERROR_FILE_NOT_FOUND = 200
    integer, parameter :: ERROR_FILE_READ = 201
    integer, parameter :: ERROR_FILE_WRITE = 202

    ! Error codes - Memory (300-399)
    integer, parameter :: ERROR_ALLOCATION_FAILED = 300
    integer, parameter :: ERROR_DEALLOCATION_FAILED = 301

contains

    !> Get human-readable error message for a given error code
    !!
    !! @param[in] ierr  Error code
    !! @return    msg   Error message string
    function get_error_message(ierr) result(msg)
        integer, intent(in) :: ierr
        character(len=256) :: msg

        select case (ierr)
        ! Success
        case (ERROR_SUCCESS)
            msg = "Success"

        ! Input validation errors (1-99)
        case (ERROR_INVALID_INPUT)
            msg = "Invalid input parameter"
        case (ERROR_OUT_OF_BOUNDS)
            msg = "Value out of valid bounds"
        case (ERROR_SIZE_MISMATCH)
            msg = "Array size mismatch"
        case (ERROR_INVALID_CONCENTRATION)
            msg = "Invalid concentration (must be 0 < c <= 100)"
        case (ERROR_NEGATIVE_VALUE)
            msg = "Value must be non-negative"
        case (ERROR_INVALID_RANGE)
            msg = "Value outside valid range"

        ! Numerical errors (100-199)
        case (ERROR_CONVERGENCE_FAILED)
            msg = "Numerical convergence failed"
        case (ERROR_SINGULAR_MATRIX)
            msg = "Matrix is singular or ill-conditioned"
        case (ERROR_LAPACK_FAILED)
            msg = "LAPACK routine failed"

        ! I/O errors (200-299)
        case (ERROR_FILE_NOT_FOUND)
            msg = "File not found"
        case (ERROR_FILE_READ)
            msg = "Error reading file"
        case (ERROR_FILE_WRITE)
            msg = "Error writing file"

        ! Memory errors (300-399)
        case (ERROR_ALLOCATION_FAILED)
            msg = "Memory allocation failed"
        case (ERROR_DEALLOCATION_FAILED)
            msg = "Memory deallocation failed"

        case default
            write(msg, '(A,I0)') "Unknown error code: ", ierr
        end select
    end function get_error_message

    !> Generic error handler that prints error messages
    !!
    !! @param[in] ierr          Error code
    !! @param[in] routine_name  Name of the routine where error occurred
    !! @param[in] extra_info    Optional additional information (default: "")
    !! @param[in] fatal         If .true., stop execution (default: .false.)
    subroutine error_handler(ierr, routine_name, extra_info, fatal)
        integer, intent(in) :: ierr
        character(len=*), intent(in) :: routine_name
        character(len=*), intent(in), optional :: extra_info
        logical, intent(in), optional :: fatal
        character(len=256) :: msg
        logical :: is_fatal

        ! Check if fatal (default is .false.)
        is_fatal = .false.
        if (present(fatal)) is_fatal = fatal

        ! If success, do nothing
        if (ierr == ERROR_SUCCESS) return

        ! Get error message
        msg = get_error_message(ierr)

        ! Print error
        write(*, '(A)') "============================================"
        write(*, '(A,I0)') "ERROR CODE: ", ierr
        write(*, '(A,A)') "ROUTINE:    ", trim(routine_name)
        write(*, '(A,A)') "MESSAGE:    ", trim(msg)
        if (present(extra_info)) then
            if (len_trim(extra_info) > 0) then
                write(*, '(A,A)') "DETAILS:    ", trim(extra_info)
            end if
        end if
        write(*, '(A)') "============================================"

        ! Stop if fatal
        if (is_fatal) then
            write(*, '(A)') "FATAL ERROR - Stopping execution"
            stop 1
        end if
    end subroutine error_handler

    !> Check if a value is within bounds [min_val, max_val]
    !!
    !! @param[in]  val      Value to check
    !! @param[in]  min_val  Minimum allowed value
    !! @param[in]  max_val  Maximum allowed value
    !! @param[out] ierr     ERROR_SUCCESS or ERROR_OUT_OF_BOUNDS
    subroutine check_bounds(val, min_val, max_val, ierr)
        integer, intent(in) :: val, min_val, max_val
        integer, intent(out) :: ierr

        if (val < min_val .or. val > max_val) then
            ierr = ERROR_OUT_OF_BOUNDS
        else
            ierr = ERROR_SUCCESS
        end if
    end subroutine check_bounds

    !> Check if a value is positive (> 0) or non-negative (>= 0)
    !!
    !! @param[in]  val           Value to check
    !! @param[in]  allow_zero    If .true., accepts 0; if .false., requires > 0
    !! @param[out] ierr          ERROR_SUCCESS or ERROR_NEGATIVE_VALUE
    subroutine check_positive(val, allow_zero, ierr)
        real(dp), intent(in) :: val
        logical, intent(in) :: allow_zero
        integer, intent(out) :: ierr

        if (allow_zero) then
            ! Check val >= 0
            if (val < 0.0_dp) then
                ierr = ERROR_NEGATIVE_VALUE
            else
                ierr = ERROR_SUCCESS
            end if
        else
            ! Check val > 0
            if (val <= 0.0_dp) then
                ierr = ERROR_NEGATIVE_VALUE
            else
                ierr = ERROR_SUCCESS
            end if
        end if
    end subroutine check_positive

    !> Check if a value is within a real range [min_val, max_val]
    !!
    !! @param[in]  val      Value to check
    !! @param[in]  min_val  Minimum allowed value
    !! @param[in]  max_val  Maximum allowed value
    !! @param[out] ierr     ERROR_SUCCESS or ERROR_INVALID_RANGE
    subroutine check_range(val, min_val, max_val, ierr)
        real(dp), intent(in) :: val, min_val, max_val
        integer, intent(out) :: ierr

        if (val < min_val .or. val > max_val) then
            ierr = ERROR_INVALID_RANGE
        else
            ierr = ERROR_SUCCESS
        end if
    end subroutine check_range

end module lsda_errors