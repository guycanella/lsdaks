!> Unit tests for lsda_errors module
program test_lsda_errors
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp
    implicit none

    real(dp), parameter :: TOL = 1.0e-9_dp

    call execute_serial_cmd_app(get_error_tests())

contains

    function get_error_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests

        tests = test_list([ &
            test("error_codes_defined", test_error_codes_defined), &
            test("error_message_success", test_error_message_success), &
            test("error_message_input_validation", test_error_message_input_validation), &
            test("error_message_numerical", test_error_message_numerical), &
            test("error_message_io", test_error_message_io), &
            test("error_message_memory", test_error_message_memory), &
            test("error_message_unknown", test_error_message_unknown), &
            test("check_bounds_valid", test_check_bounds_valid), &
            test("check_bounds_invalid", test_check_bounds_invalid), &
            test("check_positive_allow_zero", test_check_positive_allow_zero), &
            test("check_positive_strict", test_check_positive_strict), &
            test("check_range_valid", test_check_range_valid), &
            test("check_range_invalid", test_check_range_invalid) &
        ])
    end function get_error_tests

    !> Test that all error codes are properly defined
    subroutine test_error_codes_defined()
        use fortuno_serial, only: check => serial_check
        use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT, ERROR_OUT_OF_BOUNDS, &
                               ERROR_SIZE_MISMATCH, ERROR_INVALID_CONCENTRATION, &
                               ERROR_NEGATIVE_VALUE, ERROR_INVALID_RANGE, &
                               ERROR_CONVERGENCE_FAILED, ERROR_SINGULAR_MATRIX, &
                               ERROR_LAPACK_FAILED, ERROR_FILE_NOT_FOUND, &
                               ERROR_FILE_READ, ERROR_FILE_WRITE, &
                               ERROR_ALLOCATION_FAILED, ERROR_DEALLOCATION_FAILED

        ! Test error codes are in correct ranges
        call check(ERROR_SUCCESS == 0, "ERROR_SUCCESS should be 0")

        ! Input validation errors (1-99)
        call check(ERROR_INVALID_INPUT >= 1 .and. ERROR_INVALID_INPUT <= 99, &
                   "ERROR_INVALID_INPUT in range [1,99]")
        call check(ERROR_OUT_OF_BOUNDS >= 1 .and. ERROR_OUT_OF_BOUNDS <= 99, &
                   "ERROR_OUT_OF_BOUNDS in range [1,99]")
        call check(ERROR_SIZE_MISMATCH >= 1 .and. ERROR_SIZE_MISMATCH <= 99, &
                   "ERROR_SIZE_MISMATCH in range [1,99]")
        call check(ERROR_INVALID_CONCENTRATION >= 1 .and. ERROR_INVALID_CONCENTRATION <= 99, &
                   "ERROR_INVALID_CONCENTRATION in range [1,99]")
        call check(ERROR_NEGATIVE_VALUE >= 1 .and. ERROR_NEGATIVE_VALUE <= 99, &
                   "ERROR_NEGATIVE_VALUE in range [1,99]")
        call check(ERROR_INVALID_RANGE >= 1 .and. ERROR_INVALID_RANGE <= 99, &
                   "ERROR_INVALID_RANGE in range [1,99]")

        ! Numerical errors (100-199)
        call check(ERROR_CONVERGENCE_FAILED >= 100 .and. ERROR_CONVERGENCE_FAILED <= 199, &
                   "ERROR_CONVERGENCE_FAILED in range [100,199]")
        call check(ERROR_SINGULAR_MATRIX >= 100 .and. ERROR_SINGULAR_MATRIX <= 199, &
                   "ERROR_SINGULAR_MATRIX in range [100,199]")
        call check(ERROR_LAPACK_FAILED >= 100 .and. ERROR_LAPACK_FAILED <= 199, &
                   "ERROR_LAPACK_FAILED in range [100,199]")

        ! I/O errors (200-299)
        call check(ERROR_FILE_NOT_FOUND >= 200 .and. ERROR_FILE_NOT_FOUND <= 299, &
                   "ERROR_FILE_NOT_FOUND in range [200,299]")
        call check(ERROR_FILE_READ >= 200 .and. ERROR_FILE_READ <= 299, &
                   "ERROR_FILE_READ in range [200,299]")
        call check(ERROR_FILE_WRITE >= 200 .and. ERROR_FILE_WRITE <= 299, &
                   "ERROR_FILE_WRITE in range [200,299]")

        ! Memory errors (300-399)
        call check(ERROR_ALLOCATION_FAILED >= 300 .and. ERROR_ALLOCATION_FAILED <= 399, &
                   "ERROR_ALLOCATION_FAILED in range [300,399]")
        call check(ERROR_DEALLOCATION_FAILED >= 300 .and. ERROR_DEALLOCATION_FAILED <= 399, &
                   "ERROR_DEALLOCATION_FAILED in range [300,399]")
    end subroutine test_error_codes_defined

    !> Test success error message
    subroutine test_error_message_success()
        use fortuno_serial, only: check => serial_check
        use lsda_errors, only: ERROR_SUCCESS, get_error_message

        character(len=256) :: msg

        msg = get_error_message(ERROR_SUCCESS)
        call check(trim(msg) == "Success", "Success message should be correct")
    end subroutine test_error_message_success

    !> Test input validation error messages
    subroutine test_error_message_input_validation()
        use fortuno_serial, only: check => serial_check
        use lsda_errors, only: ERROR_INVALID_INPUT, ERROR_OUT_OF_BOUNDS, &
                               ERROR_SIZE_MISMATCH, ERROR_INVALID_CONCENTRATION, &
                               ERROR_NEGATIVE_VALUE, ERROR_INVALID_RANGE, &
                               get_error_message

        character(len=256) :: msg

        msg = get_error_message(ERROR_INVALID_INPUT)
        call check(len_trim(msg) > 0, "ERROR_INVALID_INPUT should have message")

        msg = get_error_message(ERROR_OUT_OF_BOUNDS)
        call check(len_trim(msg) > 0, "ERROR_OUT_OF_BOUNDS should have message")

        msg = get_error_message(ERROR_SIZE_MISMATCH)
        call check(len_trim(msg) > 0, "ERROR_SIZE_MISMATCH should have message")

        msg = get_error_message(ERROR_INVALID_CONCENTRATION)
        call check(len_trim(msg) > 0, "ERROR_INVALID_CONCENTRATION should have message")

        msg = get_error_message(ERROR_NEGATIVE_VALUE)
        call check(len_trim(msg) > 0, "ERROR_NEGATIVE_VALUE should have message")

        msg = get_error_message(ERROR_INVALID_RANGE)
        call check(len_trim(msg) > 0, "ERROR_INVALID_RANGE should have message")
    end subroutine test_error_message_input_validation

    !> Test numerical error messages
    subroutine test_error_message_numerical()
        use fortuno_serial, only: check => serial_check
        use lsda_errors, only: ERROR_CONVERGENCE_FAILED, ERROR_SINGULAR_MATRIX, &
                               ERROR_LAPACK_FAILED, get_error_message

        character(len=256) :: msg

        msg = get_error_message(ERROR_CONVERGENCE_FAILED)
        call check(len_trim(msg) > 0, "ERROR_CONVERGENCE_FAILED should have message")

        msg = get_error_message(ERROR_SINGULAR_MATRIX)
        call check(len_trim(msg) > 0, "ERROR_SINGULAR_MATRIX should have message")

        msg = get_error_message(ERROR_LAPACK_FAILED)
        call check(len_trim(msg) > 0, "ERROR_LAPACK_FAILED should have message")
    end subroutine test_error_message_numerical

    !> Test I/O error messages
    subroutine test_error_message_io()
        use fortuno_serial, only: check => serial_check
        use lsda_errors, only: ERROR_FILE_NOT_FOUND, ERROR_FILE_READ, &
                               ERROR_FILE_WRITE, get_error_message

        character(len=256) :: msg

        msg = get_error_message(ERROR_FILE_NOT_FOUND)
        call check(len_trim(msg) > 0, "ERROR_FILE_NOT_FOUND should have message")

        msg = get_error_message(ERROR_FILE_READ)
        call check(len_trim(msg) > 0, "ERROR_FILE_READ should have message")

        msg = get_error_message(ERROR_FILE_WRITE)
        call check(len_trim(msg) > 0, "ERROR_FILE_WRITE should have message")
    end subroutine test_error_message_io

    !> Test memory error messages
    subroutine test_error_message_memory()
        use fortuno_serial, only: check => serial_check
        use lsda_errors, only: ERROR_ALLOCATION_FAILED, ERROR_DEALLOCATION_FAILED, &
                               get_error_message

        character(len=256) :: msg

        msg = get_error_message(ERROR_ALLOCATION_FAILED)
        call check(len_trim(msg) > 0, "ERROR_ALLOCATION_FAILED should have message")

        msg = get_error_message(ERROR_DEALLOCATION_FAILED)
        call check(len_trim(msg) > 0, "ERROR_DEALLOCATION_FAILED should have message")
    end subroutine test_error_message_memory

    !> Test unknown error code message
    subroutine test_error_message_unknown()
        use fortuno_serial, only: check => serial_check
        use lsda_errors, only: get_error_message

        character(len=256) :: msg
        integer :: unknown_code

        unknown_code = 9999

        msg = get_error_message(unknown_code)
        call check(index(msg, "Unknown error code") > 0, &
                   "Unknown code should return appropriate message")
    end subroutine test_error_message_unknown

    !> Test check_bounds with valid value
    subroutine test_check_bounds_valid()
        use fortuno_serial, only: check => serial_check
        use lsda_errors, only: check_bounds, ERROR_SUCCESS

        integer :: ierr

        call check_bounds(5, 1, 10, ierr)
        call check(ierr == ERROR_SUCCESS, "Value 5 in [1,10] should succeed")

        call check_bounds(1, 1, 10, ierr)
        call check(ierr == ERROR_SUCCESS, "Value 1 in [1,10] should succeed")

        call check_bounds(10, 1, 10, ierr)
        call check(ierr == ERROR_SUCCESS, "Value 10 in [1,10] should succeed")
    end subroutine test_check_bounds_valid

    !> Test check_bounds with invalid value
    subroutine test_check_bounds_invalid()
        use fortuno_serial, only: check => serial_check
        use lsda_errors, only: check_bounds, ERROR_OUT_OF_BOUNDS

        integer :: ierr

        call check_bounds(0, 1, 10, ierr)
        call check(ierr == ERROR_OUT_OF_BOUNDS, "Value 0 below [1,10] should fail")

        call check_bounds(11, 1, 10, ierr)
        call check(ierr == ERROR_OUT_OF_BOUNDS, "Value 11 above [1,10] should fail")
    end subroutine test_check_bounds_invalid

    !> Test check_positive with allow_zero = .true.
    subroutine test_check_positive_allow_zero()
        use fortuno_serial, only: check => serial_check
        use lsda_errors, only: check_positive, ERROR_SUCCESS, ERROR_NEGATIVE_VALUE
        use lsda_constants, only: dp

        integer :: ierr

        ! Test val = 0.0 (should succeed with allow_zero = .true.)
        call check_positive(0.0_dp, .true., ierr)
        call check(ierr == ERROR_SUCCESS, "Zero should succeed when allow_zero = .true.")

        ! Test val = 5.0 (should succeed)
        call check_positive(5.0_dp, .true., ierr)
        call check(ierr == ERROR_SUCCESS, "Positive value should succeed")

        ! Test val = -1.0 (should fail)
        call check_positive(-1.0_dp, .true., ierr)
        call check(ierr == ERROR_NEGATIVE_VALUE, "Negative value should fail")
    end subroutine test_check_positive_allow_zero

    !> Test check_positive with allow_zero = .false.
    subroutine test_check_positive_strict()
        use fortuno_serial, only: check => serial_check
        use lsda_errors, only: check_positive, ERROR_SUCCESS, ERROR_NEGATIVE_VALUE
        use lsda_constants, only: dp

        integer :: ierr

        ! Test val = 0.0 (should fail with allow_zero = .false.)
        call check_positive(0.0_dp, .false., ierr)
        call check(ierr == ERROR_NEGATIVE_VALUE, "Zero should fail when allow_zero = .false.")

        ! Test val = 5.0 (should succeed)
        call check_positive(5.0_dp, .false., ierr)
        call check(ierr == ERROR_SUCCESS, "Positive value should succeed")

        ! Test val = -1.0 (should fail)
        call check_positive(-1.0_dp, .false., ierr)
        call check(ierr == ERROR_NEGATIVE_VALUE, "Negative value should fail")
    end subroutine test_check_positive_strict

    !> Test check_range with valid value
    subroutine test_check_range_valid()
        use fortuno_serial, only: check => serial_check
        use lsda_errors, only: check_range, ERROR_SUCCESS
        use lsda_constants, only: dp

        integer :: ierr

        call check_range(5.0_dp, 0.0_dp, 10.0_dp, ierr)
        call check(ierr == ERROR_SUCCESS, "Value 5.0 in [0,10] should succeed")

        call check_range(0.0_dp, 0.0_dp, 10.0_dp, ierr)
        call check(ierr == ERROR_SUCCESS, "Value 0.0 in [0,10] should succeed")

        call check_range(10.0_dp, 0.0_dp, 10.0_dp, ierr)
        call check(ierr == ERROR_SUCCESS, "Value 10.0 in [0,10] should succeed")
    end subroutine test_check_range_valid

    !> Test check_range with invalid value
    subroutine test_check_range_invalid()
        use fortuno_serial, only: check => serial_check
        use lsda_errors, only: check_range, ERROR_INVALID_RANGE
        use lsda_constants, only: dp

        integer :: ierr

        call check_range(-0.1_dp, 0.0_dp, 10.0_dp, ierr)
        call check(ierr == ERROR_INVALID_RANGE, "Value -0.1 below [0,10] should fail")

        call check_range(10.1_dp, 0.0_dp, 10.0_dp, ierr)
        call check(ierr == ERROR_INVALID_RANGE, "Value 10.1 above [0,10] should fail")
    end subroutine test_check_range_invalid

end program test_lsda_errors
