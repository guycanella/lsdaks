program test_input_parser
    use fortuno_serial, only: execute_serial_cmd_app
    implicit none

    call execute_serial_cmd_app(get_input_parser_tests())

contains

    function get_input_parser_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests

        tests = test_list([ &
            test("validate_valid_inputs", test_validate_valid_inputs), &
            test("validate_L_negative", test_validate_L_negative), &
            test("validate_L_zero", test_validate_L_zero), &
            test("validate_Nup_negative", test_validate_Nup_negative), &
            test("validate_Ndown_negative", test_validate_Ndown_negative), &
            test("validate_N_exceeds_L", test_validate_N_exceeds_L), &
            test("validate_U_negative", test_validate_U_negative), &
            test("validate_bc_invalid", test_validate_bc_invalid), &
            test("validate_bc_valid_open", test_validate_bc_valid_open), &
            test("validate_bc_valid_periodic", test_validate_bc_valid_periodic), &
            test("validate_bc_valid_twisted", test_validate_bc_valid_twisted), &
            test("validate_max_iter_zero", test_validate_max_iter_zero), &
            test("validate_max_iter_negative", test_validate_max_iter_negative), &
            test("validate_mixing_alpha_zero", test_validate_mixing_alpha_zero), &
            test("validate_mixing_alpha_large", test_validate_mixing_alpha_large), &
            test("validate_mixing_alpha_negative", test_validate_mixing_alpha_negative), &
            test("convert_system_params_periodic", test_convert_system_params_periodic), &
            test("convert_system_params_open", test_convert_system_params_open), &
            test("convert_system_params_twisted", test_convert_system_params_twisted), &
            test("convert_scf_params", test_convert_scf_params) &
        ])
    end function get_input_parser_tests


    !> Test validation with valid inputs
    subroutine test_validate_valid_inputs()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(input_params_t) :: inputs
        integer :: ierr

        ! Set up valid inputs
        inputs%L = 10
        inputs%Nup = 5
        inputs%Ndown = 5
        inputs%U = 4.0_dp
        inputs%bc_type = 'periodic'
        inputs%max_iter = 100
        inputs%mixing_alpha = 0.3_dp

        call validate_inputs(inputs, ierr)

        call check(ierr == ERROR_SUCCESS, "Valid inputs should pass validation")
    end subroutine test_validate_valid_inputs


    !> Test validation with negative L
    subroutine test_validate_L_negative()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_INVALID_INPUT

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = -10  ! Invalid!
        inputs%Nup = 5
        inputs%Ndown = 5
        inputs%U = 4.0_dp
        inputs%bc_type = 'periodic'
        inputs%max_iter = 100
        inputs%mixing_alpha = 0.3_dp

        call validate_inputs(inputs, ierr)

        call check(ierr == ERROR_INVALID_INPUT, "Negative L should fail")
    end subroutine test_validate_L_negative


    !> Test validation with L = 0
    subroutine test_validate_L_zero()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_INVALID_INPUT

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = 0  ! Invalid!
        inputs%Nup = 5
        inputs%Ndown = 5
        inputs%U = 4.0_dp
        inputs%bc_type = 'periodic'
        inputs%max_iter = 100
        inputs%mixing_alpha = 0.3_dp

        call validate_inputs(inputs, ierr)

        call check(ierr == ERROR_INVALID_INPUT, "L=0 should fail")
    end subroutine test_validate_L_zero


    !> Test validation with negative Nup
    subroutine test_validate_Nup_negative()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_INVALID_INPUT

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = 10
        inputs%Nup = -5  ! Invalid!
        inputs%Ndown = 5
        inputs%U = 4.0_dp
        inputs%bc_type = 'periodic'
        inputs%max_iter = 100
        inputs%mixing_alpha = 0.3_dp

        call validate_inputs(inputs, ierr)

        call check(ierr == ERROR_INVALID_INPUT, "Negative Nup should fail")
    end subroutine test_validate_Nup_negative


    !> Test validation with negative Ndown
    subroutine test_validate_Ndown_negative()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_INVALID_INPUT

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = 10
        inputs%Nup = 5
        inputs%Ndown = -3  ! Invalid!
        inputs%U = 4.0_dp
        inputs%bc_type = 'periodic'
        inputs%max_iter = 100
        inputs%mixing_alpha = 0.3_dp

        call validate_inputs(inputs, ierr)

        call check(ierr == ERROR_INVALID_INPUT, "Negative Ndown should fail")
    end subroutine test_validate_Ndown_negative


    !> Test validation with N > L
    subroutine test_validate_N_exceeds_L()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_INVALID_INPUT

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = 10
        inputs%Nup = 8
        inputs%Ndown = 5  ! Total = 13 > L
        inputs%U = 4.0_dp
        inputs%bc_type = 'periodic'
        inputs%max_iter = 100
        inputs%mixing_alpha = 0.3_dp

        call validate_inputs(inputs, ierr)

        call check(ierr == ERROR_INVALID_INPUT, "N > L should fail")
    end subroutine test_validate_N_exceeds_L


    !> Test validation with negative U
    subroutine test_validate_U_negative()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_INVALID_INPUT

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = 10
        inputs%Nup = 5
        inputs%Ndown = 5
        inputs%U = -1.0_dp  ! Invalid!
        inputs%bc_type = 'periodic'
        inputs%max_iter = 100
        inputs%mixing_alpha = 0.3_dp

        call validate_inputs(inputs, ierr)

        call check(ierr == ERROR_INVALID_INPUT, "Negative U should fail")
    end subroutine test_validate_U_negative


    !> Test validation with invalid BC type
    subroutine test_validate_bc_invalid()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_INVALID_INPUT

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = 10
        inputs%Nup = 5
        inputs%Ndown = 5
        inputs%U = 4.0_dp
        inputs%bc_type = 'invalid_bc'  ! Invalid!
        inputs%max_iter = 100
        inputs%mixing_alpha = 0.3_dp

        call validate_inputs(inputs, ierr)

        call check(ierr == ERROR_INVALID_INPUT, "Invalid BC type should fail")
    end subroutine test_validate_bc_invalid


    !> Test validation with valid BC: open
    subroutine test_validate_bc_valid_open()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = 10
        inputs%Nup = 5
        inputs%Ndown = 5
        inputs%U = 4.0_dp
        inputs%bc_type = 'open'
        inputs%max_iter = 100
        inputs%mixing_alpha = 0.3_dp

        call validate_inputs(inputs, ierr)

        call check(ierr == ERROR_SUCCESS, "BC='open' should pass")
    end subroutine test_validate_bc_valid_open


    !> Test validation with valid BC: periodic
    subroutine test_validate_bc_valid_periodic()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = 10
        inputs%Nup = 5
        inputs%Ndown = 5
        inputs%U = 4.0_dp
        inputs%bc_type = 'periodic'
        inputs%max_iter = 100
        inputs%mixing_alpha = 0.3_dp

        call validate_inputs(inputs, ierr)

        call check(ierr == ERROR_SUCCESS, "BC='periodic' should pass")
    end subroutine test_validate_bc_valid_periodic


    !> Test validation with valid BC: twisted
    subroutine test_validate_bc_valid_twisted()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = 10
        inputs%Nup = 5
        inputs%Ndown = 5
        inputs%U = 4.0_dp
        inputs%bc_type = 'twisted'
        inputs%max_iter = 100
        inputs%mixing_alpha = 0.3_dp

        call validate_inputs(inputs, ierr)

        call check(ierr == ERROR_SUCCESS, "BC='twisted' should pass")
    end subroutine test_validate_bc_valid_twisted


    !> Test validation with max_iter = 0
    subroutine test_validate_max_iter_zero()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_INVALID_INPUT

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = 10
        inputs%Nup = 5
        inputs%Ndown = 5
        inputs%U = 4.0_dp
        inputs%bc_type = 'periodic'
        inputs%max_iter = 0  ! Invalid!
        inputs%mixing_alpha = 0.3_dp

        call validate_inputs(inputs, ierr)

        call check(ierr == ERROR_INVALID_INPUT, "max_iter=0 should fail")
    end subroutine test_validate_max_iter_zero


    !> Test validation with negative max_iter
    subroutine test_validate_max_iter_negative()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_INVALID_INPUT

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = 10
        inputs%Nup = 5
        inputs%Ndown = 5
        inputs%U = 4.0_dp
        inputs%bc_type = 'periodic'
        inputs%max_iter = -10  ! Invalid!
        inputs%mixing_alpha = 0.3_dp

        call validate_inputs(inputs, ierr)

        call check(ierr == ERROR_INVALID_INPUT, "Negative max_iter should fail")
    end subroutine test_validate_max_iter_negative


    !> Test validation with mixing_alpha = 0
    subroutine test_validate_mixing_alpha_zero()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_INVALID_INPUT

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = 10
        inputs%Nup = 5
        inputs%Ndown = 5
        inputs%U = 4.0_dp
        inputs%bc_type = 'periodic'
        inputs%max_iter = 100
        inputs%mixing_alpha = 0.0_dp  ! Invalid!

        call validate_inputs(inputs, ierr)

        call check(ierr == ERROR_INVALID_INPUT, "mixing_alpha=0 should fail")
    end subroutine test_validate_mixing_alpha_zero


    !> Test validation with mixing_alpha > 1
    subroutine test_validate_mixing_alpha_large()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_INVALID_INPUT

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = 10
        inputs%Nup = 5
        inputs%Ndown = 5
        inputs%U = 4.0_dp
        inputs%bc_type = 'periodic'
        inputs%max_iter = 100
        inputs%mixing_alpha = 1.5_dp  ! Invalid!

        call validate_inputs(inputs, ierr)

        call check(ierr == ERROR_INVALID_INPUT, "mixing_alpha>1 should fail")
    end subroutine test_validate_mixing_alpha_large


    !> Test validation with negative mixing_alpha
    subroutine test_validate_mixing_alpha_negative()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_INVALID_INPUT

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = 10
        inputs%Nup = 5
        inputs%Ndown = 5
        inputs%U = 4.0_dp
        inputs%bc_type = 'periodic'
        inputs%max_iter = 100
        inputs%mixing_alpha = -0.1_dp  ! Invalid!

        call validate_inputs(inputs, ierr)

        call check(ierr == ERROR_INVALID_INPUT, "Negative mixing_alpha should fail")
    end subroutine test_validate_mixing_alpha_negative


    !> Test convert_to_system_params with periodic BC
    subroutine test_convert_system_params_periodic()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_types, only: system_params_t
        use boundary_conditions, only: BC_PERIODIC
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(input_params_t) :: inputs
        type(system_params_t) :: sys_params
        integer :: ierr

        inputs%L = 20
        inputs%Nup = 8
        inputs%Ndown = 6
        inputs%U = 5.5_dp
        inputs%bc_type = 'periodic'
        inputs%phase = 0.5_dp

        call convert_to_system_params(inputs, sys_params, ierr)

        call check(ierr == ERROR_SUCCESS, "Conversion should succeed")
        call check(sys_params%L == 20, "L should match")
        call check(sys_params%Nup == 8, "Nup should match")
        call check(sys_params%Ndown == 6, "Ndown should match")
        call check(abs(sys_params%U - 5.5_dp) < 1.0e-10_dp, "U should match")
        call check(sys_params%bc == BC_PERIODIC, "BC should be PERIODIC")
        call check(abs(sys_params%phase - 0.5_dp) < 1.0e-10_dp, "Phase should match")
    end subroutine test_convert_system_params_periodic


    !> Test convert_to_system_params with open BC
    subroutine test_convert_system_params_open()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_types, only: system_params_t
        use boundary_conditions, only: BC_OPEN
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(input_params_t) :: inputs
        type(system_params_t) :: sys_params
        integer :: ierr

        inputs%L = 15
        inputs%Nup = 5
        inputs%Ndown = 5
        inputs%U = 3.0_dp
        inputs%bc_type = 'open'
        inputs%phase = 0.0_dp

        call convert_to_system_params(inputs, sys_params, ierr)

        call check(ierr == ERROR_SUCCESS, "Conversion should succeed")
        call check(sys_params%bc == BC_OPEN, "BC should be OPEN")
    end subroutine test_convert_system_params_open


    !> Test convert_to_system_params with twisted BC
    subroutine test_convert_system_params_twisted()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_types, only: system_params_t
        use boundary_conditions, only: BC_TWISTED
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(input_params_t) :: inputs
        type(system_params_t) :: sys_params
        integer :: ierr

        inputs%L = 12
        inputs%Nup = 4
        inputs%Ndown = 4
        inputs%U = 2.0_dp
        inputs%bc_type = 'twisted'
        inputs%phase = 1.0_dp

        call convert_to_system_params(inputs, sys_params, ierr)

        call check(ierr == ERROR_SUCCESS, "Conversion should succeed")
        call check(sys_params%bc == BC_TWISTED, "BC should be TWISTED")
        call check(abs(sys_params%phase - 1.0_dp) < 1.0e-10_dp, "Phase should match")
    end subroutine test_convert_system_params_twisted


    !> Test convert_to_scf_params
    subroutine test_convert_scf_params()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use kohn_sham_cycle, only: scf_params_t
        use lsda_constants, only: dp

        type(input_params_t) :: inputs
        type(scf_params_t) :: scf_params

        inputs%max_iter = 200
        inputs%density_tol = 1.0e-8_dp
        inputs%energy_tol = 1.0e-10_dp
        inputs%mixing_alpha = 0.25_dp
        inputs%verbose = .false.
        inputs%store_history = .false.

        call convert_to_scf_params(inputs, scf_params)

        call check(scf_params%max_iter == 200, "max_iter should match")
        call check(abs(scf_params%density_tol - 1.0e-8_dp) < 1.0e-15_dp, &
                   "density_tol should match")
        call check(abs(scf_params%energy_tol - 1.0e-10_dp) < 1.0e-15_dp, &
                   "energy_tol should match")
        call check(abs(scf_params%mixing_alpha - 0.25_dp) < 1.0e-10_dp, &
                   "mixing_alpha should match")
        call check(.not. scf_params%verbose, "verbose should be false")
        call check(.not. scf_params%store_history, "store_history should be false")
    end subroutine test_convert_scf_params

end program test_input_parser
