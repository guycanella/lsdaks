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
            test("validate_N_exceeds_2L", test_validate_N_exceeds_2L), &
            test("validate_N_per_spin_exceeds_L", test_validate_N_per_spin_exceeds_L), &
            test("validate_U_negative_accepted", test_validate_U_negative_accepted), &
            test("validate_bc_invalid", test_validate_bc_invalid), &
            test("validate_bc_valid_open", test_validate_bc_valid_open), &
            test("validate_bc_valid_periodic", test_validate_bc_valid_periodic), &
            test("validate_bc_valid_twisted", test_validate_bc_valid_twisted), &
            test("validate_max_iter_zero", test_validate_max_iter_zero), &
            test("validate_max_iter_negative", test_validate_max_iter_negative), &
            test("validate_mixing_alpha_zero", test_validate_mixing_alpha_zero), &
            test("validate_mixing_alpha_large", test_validate_mixing_alpha_large), &
            test("validate_mixing_alpha_negative", test_validate_mixing_alpha_negative), &
            test("validate_xc_smoothing_width", test_validate_xc_smoothing_width), &
            test("validate_tolerances_positive", test_validate_tolerances_positive), &
            test("validate_tolerances_non_finite", test_validate_tolerances_non_finite), &
            test("validate_xc_smoothing_width_non_finite", &
                 test_validate_xc_smoothing_width_non_finite), &
            test("convert_system_params_periodic", test_convert_system_params_periodic), &
            test("convert_system_params_open", test_convert_system_params_open), &
            test("convert_system_params_twisted", test_convert_system_params_twisted), &
            test("convert_scf_params", test_convert_scf_params), &
            test("validate_phase_units_of_pi", test_validate_phase_units_of_pi), &
            test("parse_int_list_commas", test_parse_int_list_commas), &
            test("parse_int_list_blanks", test_parse_int_list_blanks), &
            test("parse_int_list_empty", test_parse_int_list_empty), &
            test("parse_int_list_separators", test_parse_int_list_separators), &
            test("parse_int_list_not_a_number", test_parse_int_list_not_a_number), &
            test("parse_int_list_out_of_range", test_parse_int_list_out_of_range), &
            test("parse_int_list_duplicates", test_parse_int_list_duplicates), &
            test("barrier_single_bounds_exact_width", test_barrier_single_bounds_exact_width), &
            test("namelist_new_keys", test_namelist_new_keys), &
            test("namelist_spring_constant", test_namelist_spring_constant), &
            test("namelist_unknown_key", test_namelist_unknown_key), &
            test("namelist_malformed_value", test_namelist_malformed_value), &
            test("namelist_missing_group", test_namelist_missing_group), &
            test("namelist_obsolete_distribution", test_namelist_obsolete_distribution), &
            test("namelist_no_trailing_newline", test_namelist_no_trailing_newline), &
            test("namelist_no_closing_slash", test_namelist_no_closing_slash), &
            test("namelist_dollar_group_form", test_namelist_dollar_group_form), &
            test("namelist_imp_positions_overflow", test_namelist_imp_positions_overflow), &
            test("namelist_group_gate_matches_reader", &
                 test_namelist_group_gate_matches_reader) &
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


    !> Test the Pauli exclusion bound on the total particle number
    !!
    !! Each site holds at most one spin-up and one spin-down electron, so the
    !! physical bound is Nup + Ndown <= 2*L, NOT Nup + Ndown <= L. Double
    !! occupancy is a legitimate configuration of the Hubbard model (it is what
    !! the U term acts on), so N > L must be accepted. This test pins both
    !! sides of the real boundary: 2*L + 1 is rejected, exactly 2*L is accepted.
    !!
    !! The total bound is enforced as a consequence of the per-spin bounds
    !! (0 <= Nup <= L, 0 <= Ndown <= L), which are pinned by
    !! test_validate_N_per_spin_exceeds_L.
    subroutine test_validate_N_exceeds_2L()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_INVALID_INPUT, ERROR_SUCCESS

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = 10
        inputs%U = 4.0_dp
        inputs%bc_type = 'periodic'
        inputs%max_iter = 100
        inputs%mixing_alpha = 0.3_dp

        ! Above the Pauli bound: Nup + Ndown = 21 > 2*L = 20
        inputs%Nup = 11
        inputs%Ndown = 10
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "N > 2L should fail (Pauli exclusion)")

        ! Exactly at the Pauli bound: completely filled band, Nup + Ndown = 2*L
        inputs%Nup = 10
        inputs%Ndown = 10
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_SUCCESS, "N = 2L should be accepted (band completely filled)")

        ! Between L and 2L: double occupancy is physical and must be accepted
        inputs%Nup = 8
        inputs%Ndown = 5
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_SUCCESS, "L < N < 2L should be accepted (double occupancy)")
    end subroutine test_validate_N_exceeds_2L


    !> Test the per-spin Pauli bound 0 <= Nup <= L and 0 <= Ndown <= L
    !!
    !! Each spin channel is a separate single-particle problem with exactly L
    !! orbitals, so no spin channel can hold more than L electrons even when the
    !! total N = Nup + Ndown stays below 2*L. Before this check existed,
    !! validate_inputs accepted e.g. Nup = 11, Ndown = 0 at L = 10 (total 11 <=
    !! 20) and the configuration was only rejected much later by
    !! compute_density_spin_real, which requires n_elec <= L. That made the
    !! public validation contract inconsistent with the solver.
    !!
    !! The lower bound is 0 (not 1): a fully polarized system Ndown = 0 is a
    !! legitimate input at this level.
    subroutine test_validate_N_per_spin_exceeds_L()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_INVALID_INPUT, ERROR_SUCCESS

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = 10
        inputs%U = 4.0_dp
        inputs%bc_type = 'periodic'
        inputs%max_iter = 100
        inputs%mixing_alpha = 0.3_dp

        ! Nup = L + 1 with a total well below 2*L: must still be rejected
        inputs%Nup = 11
        inputs%Ndown = 0
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_INVALID_INPUT, &
                   "Nup > L should fail (only L spin-up orbitals)")

        ! Ndown = L + 1 with a total well below 2*L: must still be rejected
        inputs%Nup = 0
        inputs%Ndown = 11
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_INVALID_INPUT, &
                   "Ndown > L should fail (only L spin-down orbitals)")

        ! Exactly at the per-spin bound: fully polarized filled channel
        inputs%Nup = 10
        inputs%Ndown = 0
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_SUCCESS, &
                   "Nup = L, Ndown = 0 should be accepted (fully polarized)")

        inputs%Nup = 0
        inputs%Ndown = 10
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_SUCCESS, &
                   "Ndown = L, Nup = 0 should be accepted (fully polarized)")
    end subroutine test_validate_N_per_spin_exceeds_L


    !> Test that attractive (negative) U is accepted
    !!
    !! The attractive Hubbard model (U < 0) is a supported regime of this code
    !! (it is the regime shipped in input.txt), so validate_inputs must not
    !! reject it. Only genuinely invalid parameters are rejected.
    subroutine test_validate_U_negative_accepted()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = 10
        inputs%Nup = 5
        inputs%Ndown = 5
        inputs%U = -4.0_dp  ! Attractive Hubbard interaction: valid
        inputs%bc_type = 'periodic'
        inputs%max_iter = 100
        inputs%mixing_alpha = 0.3_dp

        call validate_inputs(inputs, ierr)

        call check(ierr == ERROR_SUCCESS, "Attractive U = -4 should be accepted")
    end subroutine test_validate_U_negative_accepted


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


    !> Both convergence tolerances must be strictly positive
    !!
    !! The SCF stops when residual_V < potential_tol AND
    !! |dE| < energy_tol*max(1, |E|). Both quantities are non-negative, so a
    !! null or negative tolerance is unreachable by construction: the run would
    !! burn the whole iteration budget and be reported as a convergence failure,
    !! hiding the fact that the request itself was impossible. energy_tol in
    !! particular was NOT validated at all, even though T3 promoted it to a
    !! convergence criterion.
    !!
    !! Each boundary is probed one at a time, with everything else legal, so
    !! neither branch can be deleted without a failure here. The exact zero is
    !! included because `<= 0` and `< 0` differ precisely there.
    subroutine test_validate_tolerances_positive()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = 10
        inputs%Nup = 5
        inputs%Ndown = 5
        inputs%U = 4.0_dp
        inputs%bc_type = 'periodic'
        inputs%max_iter = 100
        inputs%mixing_alpha = 0.3_dp

        ! Baseline: the defaults are legal, so any rejection below comes from
        ! the field under test.
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_SUCCESS, &
                   "Precondition: the default tolerances must be accepted")

        inputs%potential_tol = 0.0_dp
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "potential_tol = 0 must be rejected")

        inputs%potential_tol = -1.0e-6_dp
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "negative potential_tol must be rejected")

        inputs%potential_tol = 1.0e-6_dp
        inputs%energy_tol = 0.0_dp
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "energy_tol = 0 must be rejected")

        inputs%energy_tol = -1.0e-8_dp
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "negative energy_tol must be rejected")

        ! And the smallest positive values must still pass: the check is a sign
        ! check, not a magnitude policy.
        inputs%energy_tol = tiny(1.0_dp)
        inputs%potential_tol = tiny(1.0_dp)
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_SUCCESS, &
                   "strictly positive tolerances must be accepted however small")
    end subroutine test_validate_tolerances_positive

    !> Non-finite convergence tolerances must be rejected at the input level
    !!
    !! The sign checks alone let NaN and +Infinity through:
    !!   * NaN makes its convergence comparison always false, so the criterion is
    !!     unreachable and the run only finds out after max_iter iterations;
    !!   * +Infinity makes it always true for finite residuals and energies,
    !!     which SWITCHES THAT CRITERION OFF. With potential_tol = +Inf the SCF
    !!     would declare convergence on the energy alone, i.e. exactly the false
    !!     convergence the potential-residual criterion was introduced to remove.
    !! -Infinity is already caught by `<= 0`, and is asserted here to keep that
    !! coverage explicit.
    !!
    !! Each field is probed alone with the other one legal, so deleting the
    !! finiteness guard of either tolerance makes only its own assertions fail.
    !! The special values come from ieee_value rather than from expressions like
    !! 0.0/0.0, which may be constant-folded or trap depending on the flags.
    subroutine test_validate_tolerances_non_finite()
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
                                                 ieee_positive_inf, ieee_negative_inf
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT

        type(input_params_t) :: inputs
        integer :: ierr
        real(dp) :: nan_v, pinf_v, ninf_v

        nan_v = ieee_value(1.0_dp, ieee_quiet_nan)
        pinf_v = ieee_value(1.0_dp, ieee_positive_inf)
        ninf_v = ieee_value(1.0_dp, ieee_negative_inf)

        inputs%L = 10
        inputs%Nup = 5
        inputs%Ndown = 5
        inputs%U = 4.0_dp
        inputs%bc_type = 'periodic'
        inputs%max_iter = 100
        inputs%mixing_alpha = 0.3_dp
        inputs%potential_tol = 1.0e-6_dp
        inputs%energy_tol = 1.0e-8_dp

        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_SUCCESS, &
                   "Precondition: the finite tolerances must be accepted")

        ! --- potential_tol ----------------------------------------------------
        inputs%potential_tol = nan_v
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "NaN potential_tol must be rejected")

        inputs%potential_tol = pinf_v
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "+Infinity potential_tol must be rejected")

        inputs%potential_tol = ninf_v
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "-Infinity potential_tol must be rejected")

        ! --- energy_tol -------------------------------------------------------
        inputs%potential_tol = 1.0e-6_dp
        inputs%energy_tol = nan_v
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "NaN energy_tol must be rejected")

        inputs%energy_tol = pinf_v
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "+Infinity energy_tol must be rejected")

        inputs%energy_tol = ninf_v
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "-Infinity energy_tol must be rejected")

        ! Large but finite tolerances remain legal: this is a finiteness check,
        ! not a magnitude policy.
        inputs%potential_tol = huge(1.0_dp)
        inputs%energy_tol = huge(1.0_dp)
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_SUCCESS, &
                   "finite tolerances must be accepted however large")
    end subroutine test_validate_tolerances_non_finite

    !> A NaN xc_smoothing_width must not survive the configuration path
    !!
    !! NaN passes both range comparisons and is then forwarded by app/main.f90
    !! to xc_lsda_init, where get_vxc reads `w > 0` as false and silently takes
    !! the unsmoothed branch: the user asks for smoothing and receives none.
    !! Both infinities are already rejected by the range (+Inf >= 1, -Inf < 0)
    !! and are asserted here so that any future widening of the range has to
    !! revisit the finiteness question.
    subroutine test_validate_xc_smoothing_width_non_finite()
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
                                                 ieee_positive_inf, ieee_negative_inf
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = 10
        inputs%Nup = 5
        inputs%Ndown = 5
        inputs%U = 4.0_dp
        inputs%bc_type = 'periodic'
        inputs%max_iter = 100

        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_SUCCESS, &
                   "Precondition: the default xc_smoothing_width must be accepted")

        inputs%xc_smoothing_width = ieee_value(1.0_dp, ieee_quiet_nan)
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "NaN xc_smoothing_width must be rejected")

        inputs%xc_smoothing_width = ieee_value(1.0_dp, ieee_positive_inf)
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "+Infinity xc_smoothing_width must be rejected")

        inputs%xc_smoothing_width = ieee_value(1.0_dp, ieee_negative_inf)
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "-Infinity xc_smoothing_width must be rejected")
    end subroutine test_validate_xc_smoothing_width_non_finite

    !> Test validation of xc_smoothing_width
    !!
    !! The default (0) must be accepted, since it is the C++-parity setting, a
    !! negative width is meaningless and a width >= 1 would push the smoothing
    !! window outside the physical range of n.
    subroutine test_validate_xc_smoothing_width()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT

        type(input_params_t) :: inputs
        integer :: ierr

        inputs%L = 10
        inputs%Nup = 5
        inputs%Ndown = 5
        inputs%U = 4.0_dp
        inputs%bc_type = 'periodic'
        inputs%max_iter = 100

        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_SUCCESS, "default xc_smoothing_width (0) must be accepted")

        inputs%xc_smoothing_width = 0.05_dp
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_SUCCESS, "a positive xc_smoothing_width below 1 must be accepted")

        inputs%xc_smoothing_width = -0.01_dp
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "negative xc_smoothing_width should fail")

        inputs%xc_smoothing_width = 1.0_dp
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "xc_smoothing_width >= 1 should fail")
    end subroutine test_validate_xc_smoothing_width


    !> Test convert_to_system_params with periodic BC
    subroutine test_convert_system_params_periodic()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_types, only: system_params_t
        use boundary_conditions, only: BC_PERIODIC
        use lsda_constants, only: dp, PI
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
        ! The unit conversion (units of pi -> radians) is unconditional: the
        ! phase is simply unused under periodic BC. Asserting the converted
        ! value here, rather than the raw input, keeps this test from pinning
        ! the copy-verbatim behaviour that made system_params_t%phase arrive in
        ! the wrong unit (see test_convert_system_params_twisted).
        call check(abs(sys_params%phase - 0.5_dp * PI) < 1.0e-10_dp, &
                   "Phase should be converted to radians even when unused")
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
    !!
    !! The conversion of the twist angle from units of π (the input) to radians
    !! (the unit of `system_params_t%phase`, the one the Hamiltonian builder and
    !! validate_bc_parameters work in) happens HERE, inside the converter.
    !!
    !! It used to be applied by the caller, in app/main.f90, right after this
    !! routine returned: the converter handed back a `system_params_t` whose
    !! `phase` was in the wrong unit for every one of its consumers, and the
    !! contract only closed because there happened to be a single caller. This
    !! test previously pinned that state (it asserted phase = 1.0 for an input
    !! of 1.0, i.e. no conversion at all), which is why it is updated rather
    !! than added to: the old expectation encoded the trap. Any new executable,
    !! example or test converting the inputs and feeding the result to the KS
    !! cycle would have run with θ = phase rad instead of θ = phase·π, and
    !! nothing downstream can tell 0.5 from 1.5708.
    subroutine test_convert_system_params_twisted()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_types, only: system_params_t
        use boundary_conditions, only: BC_TWISTED
        use lsda_constants, only: dp, PI
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
        call check(abs(sys_params%phase - PI) < 1.0e-10_dp, &
                   "phase = 1 (in units of pi) must come back as pi radians")

        ! A second value, so the check cannot be satisfied by a constant.
        inputs%phase = 0.5_dp
        call convert_to_system_params(inputs, sys_params, ierr)
        call check(abs(sys_params%phase - 0.5_dp * PI) < 1.0e-10_dp, &
                   "phase = 0.5 (in units of pi) must come back as pi/2 radians")
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


    !> The twist phase is an input in units of pi and its range is [0, 2)
    !!
    !! Physics: the Aharonov-Bohm phase enters the Hamiltonian as e^{i theta}
    !! with theta = phase*pi radians, so the whole physically distinct range of
    !! the twist is phase in [0, 2). Anything outside it is either redundant or
    !! rejected by validate_bc_parameters deep inside the Hamiltonian builder,
    !! where the error message can no longer name the input that caused it.
    !!
    !! Regression: before the unit convention was fixed, `phase` was consumed as
    !! if it were already in radians and no range check existed at input level,
    !! so `phase = 2.5` (= 2.5 pi) travelled all the way to the diagonalization.
    subroutine test_validate_phase_units_of_pi()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT

        type(input_params_t) :: inputs
        integer :: ierr

        inputs = input_params_t()
        inputs%L = 10
        inputs%Nup = 5
        inputs%Ndown = 5
        inputs%bc_type = 'twisted'

        inputs%phase = 0.5_dp
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_SUCCESS, "phase = 0.5 pi must be accepted")

        inputs%phase = 0.0_dp
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_SUCCESS, "phase = 0 must be accepted (reduces to periodic)")

        inputs%phase = 1.999_dp
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_SUCCESS, "phase just below 2 pi must be accepted")

        inputs%phase = 2.0_dp
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "phase = 2 (= 2 pi) must be rejected")

        inputs%phase = -0.1_dp
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "negative phase must be rejected")

        ! Under open/periodic BC the phase is unused and must not gate the run.
        inputs%bc_type = 'periodic'
        inputs%phase = 7.0_dp
        call validate_inputs(inputs, ierr)
        call check(ierr == ERROR_SUCCESS, "phase is irrelevant without twisted BC")
    end subroutine test_validate_phase_units_of_pi


    !> A comma separated list of impurity sites is parsed in order
    subroutine test_parse_int_list_commas()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_errors, only: ERROR_SUCCESS

        integer, allocatable :: values(:)
        integer :: ierr

        call parse_int_list('10, 25,40 ,55', 60, values, ierr)

        call check(ierr == ERROR_SUCCESS, "A well formed list must be accepted")
        call check(allocated(values), "values must be allocated on success")
        call check(size(values) == 4, "All four sites must be parsed")
        call check(all(values == [10, 25, 40, 55]), "Sites must keep the order given")
    end subroutine test_parse_int_list_commas


    !> Blanks alone also separate sites
    subroutine test_parse_int_list_blanks()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_errors, only: ERROR_SUCCESS

        integer, allocatable :: values(:)
        integer :: ierr

        call parse_int_list('  3   7  ', 10, values, ierr)

        call check(ierr == ERROR_SUCCESS, "Blank separated sites must be accepted")
        call check(size(values) == 2, "Two sites must be parsed")
        call check(all(values == [3, 7]), "Values must match")
    end subroutine test_parse_int_list_blanks


    !> An empty list is an error, not an empty impurity set
    !!
    !! Silently accepting it would run a clean system under the name of a
    !! disordered one.
    subroutine test_parse_int_list_empty()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_errors, only: ERROR_INVALID_INPUT

        integer, allocatable :: values(:)
        integer :: ierr

        call parse_int_list('', 10, values, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "An empty string must be rejected")
        call check(.not. allocated(values), "Nothing must be returned on error")

        call parse_int_list('     ', 10, values, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "A blank-only string must be rejected")
    end subroutine test_parse_int_list_empty


    !> Empty fields (repeated, leading or trailing commas) are errors
    !!
    !! Regression: a tolerant parser would read '10,,25' as two sites and
    !! '10,25,' as two sites as well, which is indistinguishable from the
    !! intended list - until the typo hides a site the user believed was there.
    subroutine test_parse_int_list_separators()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_errors, only: ERROR_INVALID_INPUT

        integer, allocatable :: values(:)
        integer :: ierr

        call parse_int_list('10,,25', 30, values, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "Two consecutive commas must be rejected")

        call parse_int_list(',10,25', 30, values, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "A leading comma must be rejected")

        call parse_int_list('10,25,', 30, values, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "A trailing comma must be rejected")

        call parse_int_list('10, , 25', 30, values, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "An empty field must be rejected")
    end subroutine test_parse_int_list_separators


    !> Non integer tokens are errors, never partially read numbers
    !!
    !! Regression: list-directed READ accepts '1.5' (as 1) and stops at the
    !! first bad character of '12a' (as 12), so a typo would be accepted as a
    !! different, perfectly plausible lattice site.
    subroutine test_parse_int_list_not_a_number()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_errors, only: ERROR_INVALID_INPUT

        integer, allocatable :: values(:)
        integer :: ierr

        call parse_int_list('10, abc', 30, values, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "A word must be rejected")

        call parse_int_list('12a', 30, values, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "A trailing letter must be rejected")

        call parse_int_list('1.5', 30, values, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "A real number must be rejected")

        call parse_int_list('-', 30, values, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "A lone sign must be rejected")
    end subroutine test_parse_int_list_not_a_number


    !> Sites outside [1, L] are errors
    subroutine test_parse_int_list_out_of_range()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_errors, only: ERROR_INVALID_INPUT, ERROR_SUCCESS

        integer, allocatable :: values(:)
        integer :: ierr

        call parse_int_list('0', 10, values, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "Site 0 must be rejected (1-indexed lattice)")

        call parse_int_list('11', 10, values, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "A site beyond L must be rejected")

        call parse_int_list('-3', 10, values, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "A negative site must be rejected")

        call parse_int_list('1, 10', 10, values, ierr)
        call check(ierr == ERROR_SUCCESS, "Both ends of the lattice must be accepted")
    end subroutine test_parse_int_list_out_of_range


    !> A repeated site is an error
    !!
    !! potential_impurity_multiple ADDS overlapping amplitudes, so '10, 10'
    !! would silently produce an impurity of strength 2*V0 on one site.
    subroutine test_parse_int_list_duplicates()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_errors, only: ERROR_INVALID_INPUT

        integer, allocatable :: values(:)
        integer :: ierr

        call parse_int_list('10, 25, 10', 30, values, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "A repeated site must be rejected")
    end subroutine test_parse_int_list_duplicates


    !> A single barrier covers EXACTLY `width` sites, even and odd alike
    !!
    !! Regression: the executable used to compute the barrier as
    !! [position - width/2, position + width/2] with integer division on both
    !! sides, which spans 2*(width/2)+1 sites: right for odd widths, but
    !! width+1 sites for even ones. Since the transmission through a barrier
    !! decays exponentially with its width, one extra site is a physics error.
    subroutine test_barrier_single_bounds_exact_width()
        use fortuno_serial, only: check => serial_check
        use input_parser

        integer :: i_start, i_end

        call barrier_single_bounds(50, 5, i_start, i_end)
        call check(i_end - i_start + 1 == 5, "An odd width must give exactly 5 sites")
        call check(i_start == 48 .and. i_end == 52, "An odd barrier must be centred on the site")

        call barrier_single_bounds(50, 4, i_start, i_end)
        call check(i_end - i_start + 1 == 4, "An even width must give exactly 4 sites (not 5)")
        call check(i_start == 48 .and. i_end == 51, "An even barrier keeps the centre inside")

        call barrier_single_bounds(50, 1, i_start, i_end)
        call check(i_start == 50 .and. i_end == 50, "Width 1 must be the single centre site")

        call barrier_single_bounds(10, 20, i_start, i_end)
        call check(i_end - i_start + 1 == 20, "A wide even barrier must give exactly 20 sites")
    end subroutine test_barrier_single_bounds_exact_width


    !> The namelist carries the keys the new potential/table dispatch needs
    !!
    !! Regression: quasiperiodic had no keys at all (the potential was
    !! unreachable from the executable) and the table directory was hard-wired
    !! relative to the working directory.
    subroutine test_namelist_new_keys()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        character(len=*), parameter :: fname = 'test_namelist_new_keys.txt'
        type(input_params_t) :: inputs
        integer :: io_unit, ierr

        open(newunit=io_unit, file=fname, status='replace', action='write')
        write(io_unit, '(A)') "&system"
        write(io_unit, '(A)') "  L = 8"
        write(io_unit, '(A)') "  Nup = 4"
        write(io_unit, '(A)') "  Ndown = 4"
        write(io_unit, '(A)') "  table_dir = '/tmp/lsdaks_tables'"
        write(io_unit, '(A)') "/"
        write(io_unit, '(A)') "&potential"
        write(io_unit, '(A)') "  potential_type = 'quasiperiodic'"
        write(io_unit, '(A)') "  aah_lambda = 2.5"
        write(io_unit, '(A)') "  aah_beta = 0.25"
        write(io_unit, '(A)') "  aah_phi = 1.25"
        write(io_unit, '(A)') "  imp_positions_str = '2, 5'"
        write(io_unit, '(A)') "/"
        close(io_unit)

        inputs = input_params_t()
        call read_namelist_file(fname, inputs, ierr)

        call check(ierr == ERROR_SUCCESS, "The namelist must be read")
        call check(trim(inputs%table_dir) == '/tmp/lsdaks_tables', "table_dir must be read")
        call check(abs(inputs%aah_lambda - 2.5_dp) < 1.0e-12_dp, "aah_lambda must be read")
        call check(abs(inputs%aah_beta - 0.25_dp) < 1.0e-12_dp, "aah_beta must be read")
        call check(abs(inputs%aah_phi - 1.25_dp) < 1.0e-12_dp, "aah_phi must be read")
        call check(trim(inputs%imp_positions_str) == '2, 5', "imp_positions_str must be read")

        open(newunit=io_unit, file=fname, status='old')
        close(io_unit, status='delete')
    end subroutine test_namelist_new_keys


    !> `spring_constant` must reach input_params_t from the &potential group
    !!
    !! Regression for the namelist half of Bug #5 ("Harmonic Parameter Not
    !! Passed"). The field existed in input_params_t and app/main.f90 handed it
    !! to the harmonic generator, but it was never listed in the
    !! `namelist /potential/` declaration, so
    !!
    !!     &potential
    !!       potential_type = 'harmonic'
    !!       spring_constant = 0.02
    !!     /
    !!
    !! ran with the default k = 0.001 - a trap twenty times shallower than the
    !! one requested - and, with iostat ignored, said nothing. That made the
    !! harmonic trap impossible to compare against the C++ reference: the
    !! discrepancy looked exactly like a physics bug.
    !!
    !! Without the fix this test fails twice over: `spring_constant` would stay
    !! at its default AND (now that iostat is checked) the read would abort on
    !! an unmatched key.
    subroutine test_namelist_spring_constant()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        character(len=*), parameter :: fname = 'test_namelist_spring.txt'
        type(input_params_t) :: inputs
        integer :: io_unit, ierr

        open(newunit=io_unit, file=fname, status='replace', action='write')
        write(io_unit, '(A)') "&potential"
        write(io_unit, '(A)') "  potential_type = 'harmonic'"
        write(io_unit, '(A)') "  spring_constant = 0.02"
        write(io_unit, '(A)') "/"
        close(io_unit)

        inputs = input_params_t()
        call read_namelist_file(fname, inputs, ierr)

        call check(ierr == ERROR_SUCCESS, "spring_constant must be an accepted key")
        call check(abs(inputs%spring_constant - 0.02_dp) < 1.0e-12_dp, &
                   "spring_constant must be read from the namelist, not left at 0.001")

        open(newunit=io_unit, file=fname, status='old')
        close(io_unit, status='delete')
    end subroutine test_namelist_spring_constant


    !> An unknown key must abort the read instead of being ignored
    !!
    !! Regression for the ignored `iostat`: the four namelist reads used to
    !! discard their status, so a misspelled key left its parameter silently at
    !! the default and the run completed with plausible numbers for a system
    !! nobody had asked for. Without the fix, ierr comes back ERROR_SUCCESS.
    subroutine test_namelist_unknown_key()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_errors, only: ERROR_INVALID_INPUT

        character(len=*), parameter :: fname = 'test_namelist_bad_key.txt'
        type(input_params_t) :: inputs
        integer :: io_unit, ierr

        open(newunit=io_unit, file=fname, status='replace', action='write')
        write(io_unit, '(A)') "&system"
        write(io_unit, '(A)') "  L = 12"
        write(io_unit, '(A)') "  no_such_key = 3"
        write(io_unit, '(A)') "/"
        close(io_unit)

        inputs = input_params_t()
        call read_namelist_file(fname, inputs, ierr)

        call check(ierr == ERROR_INVALID_INPUT, &
                   "An unknown namelist key must be rejected, not silently ignored")

        open(newunit=io_unit, file=fname, status='old')
        close(io_unit, status='delete')
    end subroutine test_namelist_unknown_key


    !> A malformed value must abort the read
    !!
    !! gfortran reports a malformed namelist value as iostat = -1 / "End of
    !! file", the SAME status as a group that is not in the file at all, so
    !! testing `iostat > 0` alone is not enough. The parser separates the two by
    !! scanning the file for the group header; this test is what pins that
    !! behaviour, and it fails if the check is reduced to `iostat > 0`.
    subroutine test_namelist_malformed_value()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_errors, only: ERROR_INVALID_INPUT

        character(len=*), parameter :: fname = 'test_namelist_bad_value.txt'
        type(input_params_t) :: inputs
        integer :: io_unit, ierr

        open(newunit=io_unit, file=fname, status='replace', action='write')
        write(io_unit, '(A)') "&system"
        write(io_unit, '(A)') "  L = 12"
        write(io_unit, '(A)') "  U = 1.2.3"
        write(io_unit, '(A)') "/"
        close(io_unit)

        inputs = input_params_t()
        call read_namelist_file(fname, inputs, ierr)

        call check(ierr == ERROR_INVALID_INPUT, &
                   "A malformed value must be rejected even though gfortran calls it EOF")

        open(newunit=io_unit, file=fname, status='old')
        close(io_unit, status='delete')
    end subroutine test_namelist_malformed_value


    !> A group that is simply absent is legitimate and must only warn
    !!
    !! The four groups are optional by design (input_minimal.txt has &system
    !! only). Rejecting an absent group would break every short input file, so
    !! the error handling has to stop at the broken ones.
    subroutine test_namelist_missing_group()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp, MIX_ALPHA
        use lsda_errors, only: ERROR_SUCCESS

        character(len=*), parameter :: fname = 'test_namelist_no_group.txt'
        type(input_params_t) :: inputs
        integer :: io_unit, ierr

        open(newunit=io_unit, file=fname, status='replace', action='write')
        write(io_unit, '(A)') "&system"
        write(io_unit, '(A)') "  L = 12"
        write(io_unit, '(A)') "  Nup = 6"
        write(io_unit, '(A)') "  Ndown = 6"
        write(io_unit, '(A)') "/"
        close(io_unit)

        inputs = input_params_t()
        call read_namelist_file(fname, inputs, ierr)

        call check(ierr == ERROR_SUCCESS, &
                   "An input file with only &system must still be accepted")
        call check(inputs%L == 12, "The group that IS present must be applied")
        call check(trim(inputs%potential_type) == 'uniform', &
                   "The absent groups must keep their defaults")
        call check(abs(inputs%mixing_alpha - MIX_ALPHA) < 1.0e-12_dp, &
                   "An absent &scf must leave the mixing weight at its default")

        open(newunit=io_unit, file=fname, status='old')
        close(io_unit, status='delete')
    end subroutine test_namelist_missing_group


    !> The removed key `distribution` must be rejected, not accepted
    !!
    !! `distribution` was a &potential key that nothing ever read: the disorder
    !! generator is selected by potential_type alone. It was removed, which by
    !! construction turns an old input file into one carrying an unknown key -
    !! and that must be reported (with the migration hint printed by
    !! print_namelist_hint), not absorbed.
    subroutine test_namelist_obsolete_distribution()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_errors, only: ERROR_INVALID_INPUT

        character(len=*), parameter :: fname = 'test_namelist_obsolete.txt'
        type(input_params_t) :: inputs
        integer :: io_unit, ierr

        open(newunit=io_unit, file=fname, status='replace', action='write')
        write(io_unit, '(A)') "&potential"
        write(io_unit, '(A)') "  potential_type = 'random_uniform'"
        write(io_unit, '(A)') "  distribution = 'gaussian'"
        write(io_unit, '(A)') "/"
        close(io_unit)

        inputs = input_params_t()
        call read_namelist_file(fname, inputs, ierr)

        call check(ierr == ERROR_INVALID_INPUT, &
                   "The removed key 'distribution' must be reported, not ignored")

        open(newunit=io_unit, file=fname, status='old')
        close(io_unit, status='delete')
    end subroutine test_namelist_obsolete_distribution


    !> An input file whose last line is not newline-terminated must be accepted
    !!
    !! Regression. gfortran returns iostat = -1 / "End of file" for the namelist
    !! group that closes an unterminated last line, even though every value in
    !! that group was assigned correctly. Since the parser (rightly) treats
    !! "end-of-file with the group header present" as a broken group, a
    !! perfectly valid file whose last byte is the closing '/' was rejected with
    !! the same message as a malformed value. That hit `input.txt`, the most
    !! documented way of running the code (`fpm run lsdaks` with no arguments),
    !! and four of the six files in examples/, none of which any unit test
    !! covered.
    !!
    !! The file is written in stream mode on purpose: `write(unit,'(A)')` always
    !! terminates the record, so the bug cannot be reproduced with a formatted
    !! write.
    !!
    !! The check on `output_prefix` matters: &output is the LAST group, the one
    !! whose read reaches end-of-file, so it proves the values were not just
    !! tolerated but actually applied.
    subroutine test_namelist_no_trailing_newline()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        character(len=*), parameter :: fname = 'test_namelist_no_newline.txt'
        character(len=*), parameter :: NL = new_line('a')
        character(len=*), parameter :: contents = &
            "&system" // NL // &
            "  L = 14" // NL // &
            "  Nup = 7" // NL // &
            "  Ndown = 7" // NL // &
            "  U = 2.5" // NL // &
            "/" // NL // &
            "&output" // NL // &
            "  output_prefix = 'no_newline_run'" // NL // &
            "/"
        type(input_params_t) :: inputs
        integer :: io_unit, ierr

        open(newunit=io_unit, file=fname, status='replace', action='write', &
             access='stream', form='unformatted')
        write(io_unit) contents
        close(io_unit)

        inputs = input_params_t()
        call read_namelist_file(fname, inputs, ierr)

        call check(ierr == ERROR_SUCCESS, &
                   "A file whose last line has no newline must still be accepted")
        call check(inputs%L == 14, "L must be read from a file with no final newline")
        call check(abs(inputs%U - 2.5_dp) < 1.0e-12_dp, &
                   "U must be read from a file with no final newline")
        call check(trim(inputs%output_prefix) == 'no_newline_run', &
                   "The last group, the one that ends at EOF, must be applied too")

        open(newunit=io_unit, file=fname, status='old')
        close(io_unit, status='delete')
    end subroutine test_namelist_no_trailing_newline


    !> A group with no closing '/' must still be rejected
    !!
    !! Sibling of test_namelist_no_trailing_newline: tolerating the missing final
    !! newline must not be done by tolerating end-of-file in general. This is the
    !! case that only end-of-file can reveal - the group header is in the file,
    !! the keys are well formed, and the terminator is missing - so if the
    !! end-of-file branch is ever weakened to "accept whatever was read", this
    !! test fails while the one above keeps passing.
    subroutine test_namelist_no_closing_slash()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_errors, only: ERROR_INVALID_INPUT

        character(len=*), parameter :: fname = 'test_namelist_no_slash.txt'
        type(input_params_t) :: inputs
        integer :: io_unit, ierr

        open(newunit=io_unit, file=fname, status='replace', action='write')
        write(io_unit, '(A)') "&system"
        write(io_unit, '(A)') "  L = 12"
        write(io_unit, '(A)') "  Nup = 6"
        close(io_unit)

        inputs = input_params_t()
        call read_namelist_file(fname, inputs, ierr)

        call check(ierr == ERROR_INVALID_INPUT, &
                   "A namelist group with no closing '/' must be rejected")

        open(newunit=io_unit, file=fname, status='old')
        close(io_unit, status='delete')
    end subroutine test_namelist_no_closing_slash


    !> The legacy `$group ... $end` form must be read, not silently defaulted
    !!
    !! Regression, and the worst kind: gfortran accepts the legacy header
    !! character (measured: `$system / L = 42 / $end` returns iostat = 0 with
    !! every value assigned), but `namelist_group_present` recognised only '&'.
    !! Once that function started GATING the read, a file written in the '$'
    !! form stopped being read at all: the group was reported as absent and
    !! every one of its keys went back to the default - L, Nup, Ndown, U and bc
    !! quietly replaced by 10/5/5/4.0/'periodic'. The run then completes and
    !! prints plausible numbers for a system nobody asked for, which is exactly
    !! the failure mode the iostat checking was introduced to remove.
    !!
    !! Without the fix every value check below fails (the defaults are read
    !! instead), while `ierr` still comes back ERROR_SUCCESS - hence the checks
    !! are on the values, not on the status.
    subroutine test_namelist_dollar_group_form()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        character(len=*), parameter :: fname = 'test_namelist_dollar.txt'
        type(input_params_t) :: inputs
        integer :: io_unit, ierr

        open(newunit=io_unit, file=fname, status='replace', action='write')
        write(io_unit, '(A)') "$system"
        write(io_unit, '(A)') "  L = 16"
        write(io_unit, '(A)') "  Nup = 8"
        write(io_unit, '(A)') "  Ndown = 8"
        write(io_unit, '(A)') "  U = 3.0"
        write(io_unit, '(A)') "  bc = 'open'"
        write(io_unit, '(A)') "$end"
        write(io_unit, '(A)') "$output"
        write(io_unit, '(A)') "  output_prefix = 'dollar_run'"
        write(io_unit, '(A)') "$end"
        close(io_unit)

        inputs = input_params_t()
        call read_namelist_file(fname, inputs, ierr)

        call check(ierr == ERROR_SUCCESS, "A file in the '$' form must be accepted")
        call check(inputs%L == 16, "L must be read from a '$system' group")
        call check(inputs%Nup == 8, "Nup must be read from a '$system' group")
        call check(inputs%Ndown == 8, "Ndown must be read from a '$system' group")
        call check(abs(inputs%U - 3.0_dp) < 1.0e-12_dp, &
                   "U must be read from a '$system' group, not left at the default 4.0")
        call check(trim(inputs%bc_type) == 'open', &
                   "bc must be read from a '$system' group, not left at 'periodic'")
        call check(trim(inputs%output_prefix) == 'dollar_run', &
                   "A second '$' group must be read too")

        open(newunit=io_unit, file=fname, status='old')
        close(io_unit, status='delete')
    end subroutine test_namelist_dollar_group_form


    !> A list of impurity sites too long for its field must be rejected
    !!
    !! `imp_positions_str` is a `character(len=500)` and gfortran truncates a
    !! longer namelist value SILENTLY (measured: iostat = 0, the tail simply
    !! dropped). For a list of impurity positions that is the most expensive
    !! kind of input error: the run goes on with fewer impurities than asked
    !! for and produces a perfectly plausible result for a different system.
    !! parse_int_list is deliberately strict with every other malformed list
    !! (empty field, non-integer token, out of range, repeated site), so
    !! tolerating the loss of a tail by field size was incoherent.
    !!
    !! Without the guard `read_namelist_file` returns ERROR_SUCCESS here and
    !! `imp_positions_str` holds a truncated - and syntactically valid - list.
    subroutine test_namelist_imp_positions_overflow()
        use fortuno_serial, only: check => serial_check
        use input_parser
        use lsda_errors, only: ERROR_INVALID_INPUT, ERROR_SUCCESS

        character(len=*), parameter :: fname = 'test_namelist_imp_overflow.txt'
        type(input_params_t) :: inputs
        character(len=:), allocatable :: long_list
        character(len=8) :: num
        integer :: io_unit, ierr, k

        ! 201 sites of 3 digits each ("100, 101, ...") is well past 500 chars.
        long_list = '100'
        do k = 101, 300
            write(num, '(I0)') k
            long_list = long_list // ', ' // trim(num)
        end do

        open(newunit=io_unit, file=fname, status='replace', action='write')
        write(io_unit, '(A)') "&potential"
        write(io_unit, '(A)') "  potential_type = 'impurity_multiple'"
        write(io_unit, '(A)') "  imp_positions_str = '" // trim(long_list) // "'"
        write(io_unit, '(A)') "/"
        close(io_unit)

        inputs = input_params_t()
        call read_namelist_file(fname, inputs, ierr)

        call check(ierr == ERROR_INVALID_INPUT, &
                   "A truncated imp_positions_str must be rejected, not silently shortened")

        ! A list that fits must keep working.
        open(newunit=io_unit, file=fname, status='replace', action='write')
        write(io_unit, '(A)') "&potential"
        write(io_unit, '(A)') "  imp_positions_str = '10, 25, 40'"
        write(io_unit, '(A)') "/"
        close(io_unit)

        inputs = input_params_t()
        call read_namelist_file(fname, inputs, ierr)

        call check(ierr == ERROR_SUCCESS, "A list that fits the field must be accepted")
        call check(trim(inputs%imp_positions_str) == '10, 25, 40', &
                   "A list that fits the field must be read verbatim")

        open(newunit=io_unit, file=fname, status='old')
        close(io_unit, status='delete')
    end subroutine test_namelist_imp_positions_overflow


    !> The presence gate must never disagree with what `read(nml=)` really does
    !!
    !! `namelist_group_present` decides, BEFORE the read, whether a group is in
    !! the file; a group it misses is not read at all and keeps every default in
    !! silence. That makes a divergence between the gate and the actual reader
    !! invisible in a run, so the only way to catch one is to put the two side
    !! by side on the same records - which is what this battery does. Each case
    !! asserts two things:
    !!
    !! 1. the gate returns the verdict the case documents;
    !! 2. whenever the READER assigned the sentinel value, the gate said
    !!    "present". This is the implication that matters, and it is checked
    !!    mechanically for every case rather than being reasoned about by hand:
    !!    the '$' form was accepted by the reader and missed by the gate, and a
    !!    hand-written argument that "the header is always '&'" is precisely
    !!    what let it through.
    !!
    !! The converse implication is deliberately NOT asserted: a gate false
    !! positive is harmless (the read of an absent group from an internal file
    !! returns iostat = 0 and changes nothing), and it is what a broken group
    !! legitimately looks like - present in the file, rejected by the reader.
    subroutine test_namelist_group_gate_matches_reader()
        character, parameter :: TAB = char(9)
        character, parameter :: CR = char(13)
        character(len=900) :: r(4)

        ! Tab before the header: adjustl() does not move a tab, so a gate that
        ! only looked at the first character would miss this.
        r = ''
        r(1) = TAB // '&gsys'; r(2) = '  L = 42'; r(3) = '/'
        call check_gate_vs_reader(r(1:3), 'gsys', .true., 'tab before the header')

        ! Group names are case-insensitive for the reader; the gate folds too.
        r = ''
        r(1) = '&GSYS'; r(2) = '  L = 42'; r(3) = '/'
        call check_gate_vs_reader(r(1:3), 'gsys', .true., 'uppercase header')

        r = ''
        r(1) = '&GsYs'; r(2) = '  L = 42'; r(3) = '/'
        call check_gate_vs_reader(r(1:3), 'gsys', .true., 'mixed case header')

        ! Header away from the start of the line, with and without other text.
        r = ''
        r(1) = '    &gsys L = 42 /'
        call check_gate_vs_reader(r(1:1), 'gsys', .true., 'header after blanks')

        r = ''
        r(1) = 'junk &gsys L = 42 /'
        call check_gate_vs_reader(r(1:1), 'gsys', .true., 'header after other text')

        ! Two groups on one line: both must be found.
        r = ''
        r(1) = '&gsys L = 42 / &gout M = 7 /'
        call check_gate_vs_reader(r(1:1), 'gsys', .true., 'two groups on a line (first)')
        call check_gate_vs_reader(r(1:1), 'gout', .true., 'two groups on a line (second)')

        ! A header inside a quoted VALUE is not a header: the reader ignores it
        ! and so must the gate, or an absent group would abort the run.
        r = ''
        r(1) = '&gsys'; r(2) = "  s = '&gout'"; r(3) = '  L = 42'; r(4) = '/'
        call check_gate_vs_reader(r(1:4), 'gout', .false., 'header inside single quotes')
        call check_gate_vs_reader(r(1:4), 'gsys', .true., 'the group around a quoted header')

        r = ''
        r(1) = '&gsys'; r(2) = '  s = "&gout"'; r(3) = '  L = 42'; r(4) = '/'
        call check_gate_vs_reader(r(1:4), 'gout', .false., 'header inside double quotes')

        ! A header after a comment marker is prose.
        r = ''
        r(1) = '&gsys'; r(2) = '  L = 42'; r(3) = '/'; r(4) = '! &gout M = 7 /'
        call check_gate_vs_reader(r(1:4), 'gout', .false., 'header after !')

        ! The legacy '$' form, in both its spellings. THIS is the case the gate
        ! used to miss while the reader accepted it.
        r = ''
        r(1) = '$gsys'; r(2) = '  L = 42'; r(3) = '$end'
        call check_gate_vs_reader(r(1:3), 'gsys', .true., 'legacy $group ... $end')

        r = ''
        r(1) = '$gsys'; r(2) = '  L = 42'; r(3) = '/'
        call check_gate_vs_reader(r(1:3), 'gsys', .true., 'legacy $group closed by /')

        ! '$end' is a terminator, never a group name.
        r = ''
        r(1) = '$gsys'; r(2) = '  L = 42'; r(3) = '$end'
        call check_gate_vs_reader(r(1:3), 'end', .false., '$end is not a group named end')

        ! A blank (or a tab) between the header character and the name is NOT a
        ! header for the reader - it assigns nothing - so the gate must agree.
        r = ''
        r(1) = '& gsys'; r(2) = '  L = 42'; r(3) = '/'
        call check_gate_vs_reader(r(1:3), 'gsys', .false., 'blank between & and the name')

        r = ''
        r(1) = '&' // TAB // 'gsys'; r(2) = '  L = 42'; r(3) = '/'
        call check_gate_vs_reader(r(1:3), 'gsys', .false., 'tab between & and the name')

        ! A record that still carries the CR of a CRLF file: the reader accepts
        ! it, so the gate must not be thrown off by the trailing control byte.
        r = ''
        r(1) = '&gsys' // CR; r(2) = '  L = 42' // CR; r(3) = '/' // CR
        call check_gate_vs_reader(r(1:3), 'gsys', .true., 'CRLF records')

        ! A UTF-8 BOM in front of the header: the reader skips it.
        r = ''
        r(1) = char(239) // char(187) // char(191) // '&gsys'
        r(2) = '  L = 42'; r(3) = '/'
        call check_gate_vs_reader(r(1:3), 'gsys', .true., 'UTF-8 BOM before the header')

        ! A line longer than the 512-character chunk of read_input_records.
        r = ''
        r(1) = repeat(' ', 700) // '&gsys L = 42 /'
        call check_gate_vs_reader(r(1:1), 'gsys', .true., 'header past column 512')

        ! An empty file is a single blank record (read_input_records guarantees
        ! it): no group, and the absent-group note is the right outcome.
        r = ''
        call check_gate_vs_reader(r(1:1), 'gsys', .false., 'empty file')

        r = ''
        r(1) = '! nothing here'; r(2) = '!  nor here'
        call check_gate_vs_reader(r(1:2), 'gsys', .false., 'file with comments only')

        ! Header in the last record with no body and no terminator: the group IS
        ! there and the reader fails on it (iostat = -1), which is what must be
        ! reported as a broken group instead of an absent one.
        r = ''
        r(1) = '&gsys'; r(2) = '  L = 42'; r(3) = '/'; r(4) = '&gout'
        call check_gate_vs_reader(r(1:4), 'gout', .true., 'header at the end, unterminated')
    end subroutine test_namelist_group_gate_matches_reader


    !> One case of the battery above: gate verdict vs. what the reader did
    !!
    !! The namelist groups are declared here, next to the variables they
    !! describe, so that the very same records are handed to
    !! `namelist_group_present` and to `read(records, nml=...)`.
    !!
    !! @param[in] records     Records to test (as read_input_records would yield)
    !! @param[in] group       Group to look for ('gsys', 'gout' or anything else)
    !! @param[in] expect_gate Verdict the gate must return for this case
    !! @param[in] label       Description used in the failure messages
    subroutine check_gate_vs_reader(records, group, expect_gate, label)
        use fortuno_serial, only: check => serial_check
        use input_parser, only: namelist_group_present

        character(len=*), intent(in) :: records(:)
        character(len=*), intent(in) :: group
        logical, intent(in) :: expect_gate
        character(len=*), intent(in) :: label

        integer :: L, M, io_stat
        character(len=40) :: s
        logical :: gate, assigned

        namelist /gsys/ L, s
        namelist /gout/ M

        gate = namelist_group_present(records, group)

        L = -1
        M = -1
        s = ''
        io_stat = 0
        assigned = .false.

        if (group == 'gsys') then
            read(records, nml=gsys, iostat=io_stat)
            assigned = (L == 42)
        else if (group == 'gout') then
            read(records, nml=gout, iostat=io_stat)
            assigned = (M == 7)
        end if

        call check(gate .eqv. expect_gate, &
                   "gate verdict for &" // group // " (" // label // ")")
        call check(.not. assigned .or. gate, &
                   "the reader assigned &" // group // " but the gate called it absent (" // &
                   label // ")")
    end subroutine check_gate_vs_reader

end program test_input_parser
