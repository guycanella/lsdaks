!> Main program for LSDA-Hubbard calculations
!!
!! LSDAKS: Local Spin Density Approximation - Kohn-Sham solver
!! for the 1D Hubbard model using Bethe Ansatz-based XC functionals.
program lsdaks
    use lsda_constants, only: dp, U_SMALL
    use lsda_types, only: system_params_t
    use lsda_errors, only: ERROR_SUCCESS, ERROR_FILE_NOT_FOUND, ERROR_CONVERGENCE_FAILED
    use input_parser
    use output_writer
    use kohn_sham_cycle
    use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
    use potential_factory, only: create_potential
    use potential_seed, only: resolve_random_seed
    use potential_impurity, only: potential_impurity_single, potential_impurity_random, &
                                  potential_impurity_multiple
    use boundary_conditions, only: BC_OPEN, BC_PERIODIC, BC_TWISTED
    use table_io, only: xc_table_filename
    implicit none
    
    ! Main variables
    type(input_params_t) :: inputs
    type(system_params_t) :: sys_params
    type(scf_params_t) :: scf_params
    type(scf_results_t) :: results
    type(xc_lsda_t) :: xc_func

    real(dp), allocatable :: V_ext(:)
    real(dp), allocatable :: pot_params(:)
    real(dp), allocatable :: imp_amplitudes(:)
    integer, allocatable :: imp_positions(:)
    character(len=256) :: table_file
    !> Directory where the XC table was actually found (or the first candidate
    !! that was tried, when no table exists).
    character(len=256) :: resolved_table_dir
    integer :: ierr, seed
    integer :: i_start, i_end
    !> Outcome of the SCF cycle, kept apart from `ierr`.
    !!
    !! `ierr` is reused (and overwritten) by write_results and by the cleanup
    !! helper, so the SCF status has to survive in its own variable: without it
    !! a run that ended in ERROR_CONVERGENCE_FAILED printed "NOT CONVERGED",
    !! wrote its files and then exited with status 0, telling every script and
    !! pipeline that a scientifically invalid result was a success.
    integer :: scf_status
    !> Outcome of write_results, kept apart from `ierr` for the same reason as
    !! `scf_status`: `ierr` is overwritten by the cleanup helper before the exit
    !! status is decided.
    integer :: write_status
    logical :: table_exists

    ! Timing variables
    real :: start_time, end_time, elapsed_time
    integer :: date_values(8)
    character(len=10) :: date_str, time_str
    
    ! Start timing
    call cpu_time(start_time)
    call date_and_time(date_str, time_str, values=date_values)

    call print_banner(date_str, time_str)

    call parse_inputs(inputs, ierr)
    if (ierr == -1) then
        ! Help was requested
        stop 0
    else if (ierr /= ERROR_SUCCESS) then
        stop 1
    end if
    
    call validate_inputs(inputs, ierr)
    if (ierr /= ERROR_SUCCESS) then
        print *, "ERROR: Input validation failed"
        stop 1
    end if
    
    ! convert_to_system_params is the SINGLE point where the twist angle changes
    ! unit: `inputs%phase` is in units of π (the user-facing convention,
    ! identical to the C++ reference) and `sys_params%phase` comes back in
    ! radians, the unit apply_boundary_conditions_complex and
    ! validate_bc_parameters expect. The multiplication used to be applied here,
    ! by the caller, which meant the converter returned a system_params_t in the
    ! wrong unit for every one of its consumers and only this one call site
    ! repaired it.
    call convert_to_system_params(inputs, sys_params, ierr)
    if (ierr /= ERROR_SUCCESS) then
        print *, "ERROR: Failed to convert system parameters"
        stop 1
    end if

    call convert_to_scf_params(inputs, scf_params)
    call print_configuration(inputs, sys_params, scf_params)

    if (abs(sys_params%U) < U_SMALL) then
        ! The non-interacting gas has e_xc = V_xc = 0 exactly, so no table is
        ! needed and table lookup must not reject a valid U = 0 calculation.
        print '(A)', ""
        print '(A)', "Using exact non-interacting XC functional (U = 0)"
        call xc_lsda_init(xc_func, ierr=ierr, smoothing_width=inputs%xc_smoothing_width, &
                          u_signed=sys_params%U)
    else
        ! Check XC table using |U| (tables are symmetric)
        call check_xc_table(abs(sys_params%U), inputs%table_dir, table_file, &
                            resolved_table_dir, table_exists)

        if (.not. table_exists) then
            print '(A)', ""
            print '(A)', "=========================================="
            print '(A)', "ERROR: XC table not found!"
            print '(A)', "=========================================="
            print '(A,A)', "  Requested |U| = ", trim(real_str(abs(sys_params%U), 2))
            print '(A)', ""
            print '(A)', "Please generate the XC table first using:"
            print '(A,A)', "  fpm run generate_xc_table -- --U ", trim(real_str(abs(sys_params%U), 2))
            print '(A)', ""
            print '(A,A)', "  Looked for: ", trim(table_file)
            print '(A)', ""
            call print_available_tables(resolved_table_dir)
            stop 1
        end if
    
        print '(A)', ""
        print '(A,A)', "Loading XC table: ", trim(table_file)
        if (sys_params%U < 0.0_dp) then
            print '(A)', "  Note: Using table for |U| (attractive interaction)"
            print '(A)', "        The sign of U enters via the Shiba transformation in the XC functional"
        end if

        ! The signed U must be handed over explicitly: the table file only
        ! carries |U|, so without it the attractive run would silently use the
        ! repulsive functional (the Shiba transformation would never trigger).
        call xc_lsda_init(xc_func, table_file, ierr, smoothing_width=inputs%xc_smoothing_width, &
                          u_signed=sys_params%U)
    end if
    if (ierr /= ERROR_SUCCESS) then
        print *, "ERROR: Failed to initialize XC functional"
        stop 1
    end if

    print '(A)', "  ✓ XC functional initialized"
    if (inputs%xc_smoothing_width > 0.0_dp .and. abs(sys_params%U) >= U_SMALL) then
        print '(A,F0.4)', "  Note: V_xc discontinuity at n = 1 linearly smoothed over half-width w = ", &
                          inputs%xc_smoothing_width
        print '(A)', "        (this departs from the C++ reference, which keeps the jump)"
        ! The smoothed V_xc enters compute_total_energy only through the
        ! diagonalised V_eff; the E_xc term is the UNSMOOTHED get_exc, so the
        ! two stop being a derivative pair as soon as w > 0.
        print '(A)', "        E_xc is NOT smoothed: with w > 0, V_xc is not the functional derivative"
        print '(A)', "        of the E_xc used in the total energy; the reported energy is not variational."
    end if

    allocate(V_ext(sys_params%L))

    ! Create external potential - special handling for impurities
    select case (trim(inputs%potential_type))
    
    case ('impurity')
        ! Random impurities with concentration.
        !
        ! The seed is resolved HERE, not inside the generator: with
        ! pot_seed < 0 the generator would call random_seed() and never say
        ! what it drew, leaving the realisation - and therefore the whole
        ! result - unrecoverable. resolve_random_seed draws it explicitly so
        ! that the same integer can be printed, written to the provenance of
        ! every output file and fed back as pot_seed to replay the run.
        call resolve_random_seed(inputs%pot_seed, seed)
        inputs%effective_seed = seed
        call potential_impurity_random(inputs%V0, inputs%concentration, sys_params%L, &
                                       seed, V_ext, imp_positions, ierr)
        
        if (ierr /= ERROR_SUCCESS) then
            print *, "ERROR: Failed to create random impurity potential"
            call xc_lsda_destroy(xc_func)
            deallocate(V_ext)
            stop 1
        end if
        
        print '(A,A)', "  ✓ External potential created: random impurities"
        print '(A,F0.1,A)', "    Concentration: ", inputs%concentration, "%"
        print '(A,I0,A)', "    Number of impurities: ", size(imp_positions), " sites"
        print '(A,F0.4)', "    Impurity strength: V0 = ", inputs%V0
        call print_effective_seed(inputs%pot_seed, seed)

        deallocate(imp_positions)
    
    case ('impurity_single')
        ! Single impurity at specified position
        call potential_impurity_single(inputs%V0, nint(inputs%pot_center), &
                                       sys_params%L, V_ext, ierr)
        
        if (ierr /= ERROR_SUCCESS) then
            print *, "ERROR: Failed to create single impurity potential"
            call xc_lsda_destroy(xc_func)
            deallocate(V_ext)
            stop 1
        end if
        
        print '(A,A)', "  ✓ External potential created: single impurity"
        print '(A,I0)', "    Position: site ", nint(inputs%pot_center)
        print '(A,F0.4)', "    Strength: V0 = ", inputs%V0

    case ('impurity_multiple')
        ! Several impurities of equal amplitude V0 at the sites listed in
        ! imp_positions_str. Called directly (not through the factory) for the
        ! same reason as the two cases above: the number of impurities is not
        ! expressible in the fixed-size parameter array the factory takes.
        call parse_int_list(inputs%imp_positions_str, sys_params%L, imp_positions, ierr)

        if (ierr /= ERROR_SUCCESS) then
            print *, "ERROR: invalid imp_positions_str for potential_type = 'impurity_multiple'"
            call xc_lsda_destroy(xc_func)
            deallocate(V_ext)
            stop 1
        end if

        allocate(imp_amplitudes(size(imp_positions)))
        imp_amplitudes = inputs%V0

        call potential_impurity_multiple(imp_amplitudes, imp_positions, sys_params%L, V_ext, ierr)

        if (ierr /= ERROR_SUCCESS) then
            print *, "ERROR: Failed to create multiple impurity potential"
            deallocate(imp_amplitudes)
            deallocate(imp_positions)
            call xc_lsda_destroy(xc_func)
            deallocate(V_ext)
            stop 1
        end if

        print '(A)', "  ✓ External potential created: multiple impurities"
        print '(A,I0)', "    Number of impurities: ", size(imp_positions)
        print '(A,A)', "    Sites: ", trim(inputs%imp_positions_str)
        print '(A,A)', "    Strength: V0 = ", trim(real_str(inputs%V0, 4))

        deallocate(imp_amplitudes)
        deallocate(imp_positions)

    case default
        ! Use factory for other potential types
        select case (trim(inputs%potential_type))
        case ("uniform")
            allocate(pot_params(1))
            pot_params(1) = inputs%V0
        case ("harmonic")
            allocate(pot_params(1))
            pot_params(1) = inputs%spring_constant
        case ("random_uniform", "random_gaussian")
            allocate(pot_params(1))
            pot_params(1) = inputs%disorder_strength
        case ("barrier_single")
            allocate(pot_params(3))
            pot_params(1) = inputs%V0
            ! Exactly `width` sites, for even and odd widths alike. The old
            ! expression `position +/- width/2` used integer division on both
            ! sides and covered width+1 sites whenever width was even.
            call barrier_single_bounds(inputs%position, inputs%width, i_start, i_end)
            pot_params(2) = real(i_start, dp)
            pot_params(3) = real(i_end, dp)
        case ("barrier_double")
            allocate(pot_params(4))
            ! Matches C++ double_barrier(Na, Vb, Lb, Vwell, Lwell, v_ext)
            pot_params(1) = inputs%V0           ! V_bar: barrier height
            pot_params(2) = inputs%barrier_width ! L_bar: barrier width
            pot_params(3) = inputs%well_depth    ! V_well: well potential (typically negative)
            pot_params(4) = inputs%well_width    ! L_well: well width
        case ("quasiperiodic")
            ! Aubry-André-Harper: V(i) = lambda*cos(2*pi*beta*i + phi).
            ! Without this case the type fell into the default below, which
            ! passes a single parameter, and create_potential rejected it with
            ! ERROR_INVALID_INPUT: the potential was unreachable from the
            ! executable even though it is implemented and advertised.
            allocate(pot_params(3))
            pot_params(1) = inputs%aah_lambda
            pot_params(2) = inputs%aah_beta
            pot_params(3) = inputs%aah_phi
        case default
            allocate(pot_params(1))
            pot_params(1) = inputs%V0
        end select

        ! Same explicit resolution as in the 'impurity' branch above: the
        ! factory forwards the seed to the random generators, which would
        ! otherwise draw an unrecorded one whenever pot_seed < 0. Resolving it
        ! for every type (not only the stochastic ones) keeps this single line
        ! honest: `seed` is always the integer that was actually used.
        call resolve_random_seed(inputs%pot_seed, seed)
        inputs%effective_seed = seed

        call create_potential(inputs%potential_type, pot_params, sys_params%L, seed, V_ext, ierr)
        deallocate(pot_params)

        if (ierr /= ERROR_SUCCESS) then
            print *, "ERROR: Failed to create external potential of type '", &
                     trim(inputs%potential_type), "'"
            print '(A)', "  Supported potential_type values:"
            print '(A)', "    uniform, harmonic, quasiperiodic,"
            print '(A)', "    impurity, impurity_single, impurity_multiple,"
            print '(A)', "    random_uniform, random_gaussian,"
            print '(A)', "    barrier_single, barrier_double"
            call xc_lsda_destroy(xc_func)
            deallocate(V_ext)
            stop 1
        end if

        print '(A,A)', "  ✓ External potential created: ", trim(inputs%potential_type)
        if (trim(inputs%potential_type) == "harmonic") then
            print '(A,A)', "    spring_constant (k): ", trim(real_str(inputs%spring_constant, 6))
        end if
        if (trim(inputs%potential_type) == "random_uniform" .or. &
            trim(inputs%potential_type) == "random_gaussian") then
            print '(A,A)', "    disorder_strength: ", trim(real_str(inputs%disorder_strength, 6))
            call print_effective_seed(inputs%pot_seed, seed)
        end if
        if (trim(inputs%potential_type) == "quasiperiodic") then
            print '(A,A)', "    lambda: ", trim(real_str(inputs%aah_lambda, 4))
            print '(A,A)', "    beta:   ", trim(real_str(inputs%aah_beta, 6))
            print '(A,A)', "    phi:    ", trim(real_str(inputs%aah_phi, 4))
        end if
    end select
    
    print '(A)', ""
    
    ! NOTE: `results` is deliberately NOT initialized here. run_kohn_sham_scf_*
    ! declares it intent(out) and owns its initialization end to end (reset on
    ! entry, history allocated internally and only when store_history is set).
    ! Calling init_scf_results here would allocate a convergence history that the
    ! very next call discards.

    print '(A)', "=========================================="
    print '(A)', "  Starting Kohn-Sham SCF Cycle"
    print '(A)', "=========================================="
    print '(A)', ""
    
    if (sys_params%bc == BC_TWISTED) then
        print '(A)', "Using complex Hamiltonian (twisted BC)"
        call run_kohn_sham_scf_complex(sys_params, scf_params, V_ext, xc_func, results, ierr)
    else
        print '(A)', "Using real Hamiltonian"
        call run_kohn_sham_scf_real(sys_params, scf_params, V_ext, xc_func, results, ierr)
    end if
    
    print '(A)', ""

    ! Preserve the SCF outcome before `ierr` is reused by the writers below.
    scf_status = ierr

    if (scf_status /= ERROR_SUCCESS) then
        if (scf_status == ERROR_CONVERGENCE_FAILED) then
            print '(A)', "=========================================="
            print '(A)', "WARNING: SCF did not converge!"
            print '(A)', "=========================================="
            print '(A,I0)', "  Iterations performed: ", results%n_iterations
            print '(A,ES12.4)', "  Final |ΔV| residual:  ", results%final_potential_residual
            print '(A,ES12.4)', "  Final density error:  ", results%final_density_error
            print '(A,F16.8)', "  Final energy:         ", results%final_energy
            print '(A)', ""
            print '(A)', "The potential is NOT self-consistent; results are NOT converged."
            print '(A)', "Results may be unreliable."
            print '(A)', "Consider:"
            print '(A)', "  - Increasing max_iter"
            print '(A)', "  - Adjusting mixing_alpha"
            print '(A)', "  - Checking system parameters"
            print '(A)', ""
        else
            print *, "ERROR: SCF calculation failed with error code:", scf_status
            call cleanup_and_exit(xc_func, V_ext, results)
            stop 1
        end if
    end if
    
    print '(A)', "=========================================="
    print '(A)', "  Writing Results"
    print '(A)', "=========================================="
    print '(A)', ""
    
    call write_results(results, sys_params, inputs, ierr)
    ! Preserve the writing outcome too: `ierr` is reused by cleanup_and_exit
    ! below, and a run whose results never reached the disk must not report
    ! success. A valid calculation whose record was lost is, for every script
    ! and pipeline downstream, indistinguishable from a calculation that was
    ! never performed - unless the exit status says so.
    write_status = ierr
    if (write_status /= ERROR_SUCCESS) then
        print '(A)', "ERROR: some output files could not be written."
    end if

    call cleanup_and_exit(xc_func, V_ext, results)

    ! End timing
    call cpu_time(end_time)
    elapsed_time = end_time - start_time

    print '(A)', ""
    print '(A)', "=========================================="
    print '(A)', "  Calculation Complete!"
    print '(A)', "=========================================="
    print '(A)', ""
    print '(A,F12.3,A)', "Elapsed CPU Time: ", elapsed_time, " seconds"
    print '(A)', ""

    ! A non-self-consistent run is a FAILED run, however complete its output
    ! files look. The results were written on purpose (they are the best
    ! available estimate and the user must be able to inspect them) and the
    ! resources were released above, so the only thing left is to tell the
    ! caller the truth through the exit status.
    if (scf_status == ERROR_CONVERGENCE_FAILED) then
        print '(A)', "Exiting with status 1: the SCF cycle did NOT converge."
        print '(A)', ""
        stop 1
    end if

    ! A converged calculation whose output could not be written is also a failed
    ! run: the numbers are gone and only the exit status can say so.
    if (write_status /= ERROR_SUCCESS) then
        print '(A)', "Exiting with status 1: the results could NOT be written to disk."
        print '(A)', ""
        stop 1
    end if

contains

    !> Decimal text of a real, unpadded but with a guaranteed leading zero
    !!
    !! The `F0.d` edit descriptor produces the value with no padding, which is
    !! what these messages want, but gfortran writes it WITHOUT the integer zero
    !! when the magnitude is below one: `0.001` comes out as ".001000". The
    !! project prints "0.2000", not ".2000", and several of the values reported
    !! here are routinely below one (the harmonic k, the AAH beta), so the zero
    !! is restored explicitly instead of padding the field.
    !!
    !! @param[in] x        Value to format
    !! @param[in] decimals Number of digits after the decimal point
    !! @return    Left-justified text of `x`, trailing-blank padded
    function real_str(x, decimals) result(str)
        real(dp), intent(in) :: x
        integer, intent(in) :: decimals
        character(len=40) :: str
        character(len=16) :: fmt

        write(fmt, '(A,I0,A)') '(F0.', decimals, ')'
        write(str, fmt) x

        if (str(1:1) == '.') then
            str = '0' // trim(str)
        else if (len_trim(str) >= 2) then
            if (str(1:2) == '-.') str = '-0' // trim(str(2:))
        end if
    end function real_str

    !> Report the seed that actually produced the disorder realisation
    !!
    !! Printed for the stochastic potentials only. When the user asked for a
    !! drawn seed (`pot_seed < 0`) the message also states how to replay the
    !! run, because that is the whole point of capturing the value: the same
    !! integer appears in the provenance header of every output file.
    !!
    !! @param[in] requested Seed as given in the input (< 0 means "draw one")
    !! @param[in] effective Seed actually handed to the generator
    subroutine print_effective_seed(requested, effective)
        integer, intent(in) :: requested, effective

        print '(A,I0)', "    Random seed (effective): ", effective
        if (requested < 0) then
            print '(A,I0,A)', "      (drawn from the system clock; set pot_seed = ", &
                              effective, " to reproduce this realisation)"
        end if
    end subroutine print_effective_seed

    subroutine print_banner(date_str, time_str)
        character(len=*), intent(in) :: date_str, time_str
        character(len=50) :: formatted_date, formatted_time

        ! Format date: YYYYMMDD -> YYYY-MM-DD
        formatted_date = date_str(1:4) // '-' // date_str(5:6) // '-' // date_str(7:8)

        ! Format time: HHMMSS.sss -> HH:MM:SS
        formatted_time = time_str(1:2) // ':' // time_str(3:4) // ':' // time_str(5:6)

        print '(A)', ""
        print '(A)', "=========================================="
        print '(A)', "       LSDAKS - Hubbard LSDA Solver"
        print '(A)', "=========================================="
        print '(A)', ""
        print '(A)', "1D Hubbard Model with LSDA-DFT"
        print '(A)', "Using Bethe Ansatz XC Functionals"
        print '(A)', ""
        print '(A)', "Author: Guilherme Canella"
        print '(A)', "Version: 0.1.0"
        print '(A)', ""
        print '(A,A)', "Run Date: ", trim(formatted_date)
        print '(A,A)', "Run Time: ", trim(formatted_time)
        print '(A)', ""
        print '(A)', "=========================================="
        print '(A)', ""
    end subroutine print_banner

    subroutine print_configuration(inputs, sys_params, scf_params)
        type(input_params_t), intent(in) :: inputs
        type(system_params_t), intent(in) :: sys_params
        type(scf_params_t), intent(in) :: scf_params
        
        print '(A)', "Configuration:"
        print '(A)', "----------------------------------------"
        print '(A)', "System:"
        print '(A,I0)', "  L (sites):        ", sys_params%L
        print '(A,I0)', "  N_up:             ", sys_params%Nup
        print '(A,I0)', "  N_down:           ", sys_params%Ndown
        print '(A,F6.4)', "  Filling:          ", real(sys_params%Nup + sys_params%Ndown, dp) / real(sys_params%L, dp)
        print '(A,F0.4)', "  U:                ", sys_params%U
        if (sys_params%U < 0.0_dp) then
            print '(A)', "    (attractive interaction - pairing)"
        end if
        ! Convert BC to full name
        select case (trim(inputs%bc_type))
        case ('open')
            print '(A)', "  BC:               Open Boundary Condition"
        case ('periodic')
            print '(A)', "  BC:               Periodic Boundary Condition"
        case ('twisted')
            print '(A)', "  BC:               Twisted Boundary Condition"
        case default
            print '(A,A)', "  BC:               ", trim(inputs%bc_type)
        end select
        if (sys_params%bc == BC_TWISTED) then
            ! Printed in BOTH units on purpose: the input is in units of π (like
            ! the C++ reference, which also reports phase/pi) while the
            ! Hamiltonian uses radians, and the old line printed the radian
            ! field with a "π" suffix.
            ! F6.4 / F8.6 rather than F0.4 / F0.6: the project prints "0.2000",
            ! not ".2000", and the phase is the one configuration value that is
            ! almost always below 1 (it is bounded by 2), so the edit descriptor
            ! with no leading zero would drop it on essentially every twisted
            ! run. The widths are enough for the validated range [0, 2).
            print '(A,F6.4,A,F8.6,A)', "  Phase:            ", inputs%phase, " π  (= ", &
                                       sys_params%phase, " rad)"
        end if
        print '(A,A)', "  XC table dir:     ", trim(inputs%table_dir)
        print '(A)', ""
        
        print '(A)', "Potential:"
        print '(A,A)', "  Type:             ", trim(inputs%potential_type)
        if (abs(inputs%V0) > 1.0e-10_dp) then
            print '(A,F0.4)', "  Strength (V0):    ", inputs%V0
        end if
        if (trim(inputs%potential_type) == 'impurity') then
            print '(A,F0.1,A)', "  Concentration:    ", inputs%concentration, "%"
        else if (trim(inputs%potential_type) == 'impurity_single') then
            print '(A,I0)', "  Position:         ", nint(inputs%pot_center)
        end if
        print '(A)', ""
        
        print '(A)', "SCF Parameters:"
        print '(A,I0)', "  Max iterations:   ", scf_params%max_iter
        print '(A,ES10.2)', "  Potential tol:    ", scf_params%potential_tol
        print '(A,ES10.2)', "  Energy tol:       ", scf_params%energy_tol
        ! density_tol is NOT a convergence criterion any more (T3): the SCF stops
        ! on the potential residual plus energy stability. It is printed apart and
        ! explicitly labelled so nobody expects to tighten/loosen convergence with it.
        print '(A,ES10.2)', "  Density tol:      ", scf_params%density_tol
        print '(A)', "    (diagnostic only - not a convergence criterion)"
        if (scf_params%use_adaptive_mixing) then
            print '(A,F5.3)', "  Mixing alpha:     ", scf_params%mixing_alpha
            print '(A)', "    (initial value; adaptive mixing retunes it)"
        else
            print '(A,F5.3)', "  Mixing alpha:     ", scf_params%mixing_alpha
        end if
        print '(A,ES10.2)', "  XC smoothing w:   ", inputs%xc_smoothing_width
        if (inputs%xc_smoothing_width > 0.0_dp) then
            print '(A)', "    (V_xc discontinuity at n = 1 smoothed - NOT C++ parity)"
        end if
        print '(A,L1)', "  Verbose:          ", scf_params%verbose
        print '(A)', "----------------------------------------"
    end subroutine print_configuration

    !> Locate the XC table for a given U (uses |U|, the tables are symmetric)
    !!
    !! The table directory is no longer hard-wired relative to the working
    !! directory: the executable must be usable from anywhere. Two candidates
    !! are tried, in order:
    !!
    !! 1. `table_dir` (the `table_dir` key of the `&system` namelist, default
    !!    `data/tables/fortran_native`);
    !! 2. the directory named by the environment variable `LSDAKS_TABLE_DIR`,
    !!    when it is set and non-blank.
    !!
    !! @param[in]  U            Hubbard interaction (the sign is ignored)
    !! @param[in]  table_dir    Preferred table directory
    !! @param[out] table_file   Path of the table (existing one if found, else
    !!                          the path that was tried first)
    !! @param[out] resolved_dir Directory the returned path belongs to
    !! @param[out] exists       .true. if the table was found
    subroutine check_xc_table(U, table_dir, table_file, resolved_dir, exists)
        real(dp), intent(in) :: U
        character(len=*), intent(in) :: table_dir
        character(len=*), intent(out) :: table_file
        character(len=*), intent(out) :: resolved_dir
        logical, intent(out) :: exists

        character(len=256) :: env_dir, candidate
        integer :: env_len, env_status

        ! Candidate 1: the configured directory.
        resolved_dir = adjustl(table_dir)
        call xc_table_filename(resolved_dir, U, table_file)
        inquire(file=table_file, exist=exists)
        if (exists) return

        ! Candidate 2: the environment variable.
        call get_environment_variable('LSDAKS_TABLE_DIR', env_dir, env_len, env_status)
        if (env_status == 0 .and. env_len > 0) then
            candidate = adjustl(env_dir)
            call xc_table_filename(candidate, U, table_file)
            inquire(file=table_file, exist=exists)
            if (exists) then
                resolved_dir = candidate
                return
            end if
        end if

        ! Nothing found: report the first candidate, which is the one the user
        ! configured and therefore the one worth naming in the error message.
        resolved_dir = adjustl(table_dir)
        call xc_table_filename(resolved_dir, U, table_file)
        exists = .false.
    end subroutine check_xc_table


    !> Print the XC tables that are ACTUALLY available in a directory
    !!
    !! The previous version probed only the integers 1..20 and the values
    !! x.10, so every table whose name ends in another hundredth - all the
    !! x.90 tables shipped in data/tables/fortran_native, for instance - was
    !! reported as missing. The scan now covers the whole grid of two-decimal
    !! values in (0, U_SCAN_MAX], which is the only naming the writer can
    !! produce (the file name is formatted with F0.2).
    !!
    !! A real directory listing (execute_command_line + `ls`) was rejected on
    !! purpose: it would need a shell, a writable temporary file and a parser
    !! for its output, all on the error path of a failed run. Probing the
    !! finite set of representable names is pure Fortran and deterministic.
    !!
    !! @param[in] table_dir Directory to scan
    subroutine print_available_tables(table_dir)
        character(len=*), intent(in) :: table_dir
        !> Largest |U| probed. Tables above this are not listed (they are still
        !! usable: check_xc_table looks the requested one up directly).
        integer, parameter :: U_SCAN_MAX = 40
        !> Probing step, in hundredths: the file names carry two decimals.
        integer, parameter :: N_PROBES = U_SCAN_MAX * 100
        character(len=256) :: filename
        integer :: i, n_found
        real(dp) :: U_value
        real(dp), allocatable :: found_U(:)
        logical :: exists

        inquire(file=trim(table_dir), exist=exists)
        if (.not. exists) then
            print '(A,A)', "Table directory not found: ", trim(table_dir)
            print '(A)', ""
            print '(A)', "Please create the directory and generate tables using:"
            print '(A,A)', "  mkdir -p ", trim(table_dir)
            print '(A)', "  fpm run generate_xc_table -- --U <value>"
            print '(A)', ""
            print '(A)', "You can also point the code at an existing directory with the"
            print '(A)', "table_dir key of the &system namelist or with LSDAKS_TABLE_DIR."
            return
        end if

        print '(A)', "Scanning for available XC tables in:"
        print '(A,A)', "  ", trim(table_dir)
        print '(A)', ""

        allocate(found_U(N_PROBES))
        n_found = 0

        do i = 1, N_PROBES
            U_value = real(i, dp) / 100.0_dp
            call xc_table_filename(table_dir, U_value, filename)
            inquire(file=filename, exist=exists)
            if (exists) then
                n_found = n_found + 1
                found_U(n_found) = U_value
            end if
        end do

        if (n_found == 0) then
            print '(A)', "No XC tables found!"
            print '(A)', ""
            print '(A)', "Generate tables using:"
            print '(A)', "  fpm run generate_xc_table -- --U <value>"
            print '(A)', ""
        else
            print '(A,I0,A)', "Found ", n_found, " available table(s):"
            print '(A)', ""
            
            do i = 1, n_found
                if (mod(i-1, 5) == 0 .and. i > 1) print '(A)', ""
                write(*, '(A,F6.2)', advance='no') "  |U| = ", found_U(i)
            end do
            
            print '(A)', ""
            print '(A)', ""
            print '(A)', "Note: Tables work for both positive and negative U"
            print '(A)', "      (U=4 and U=-4 use the same table)"
            print '(A)', ""
        end if
        
        deallocate(found_U)
        
    end subroutine print_available_tables

    !> Cleanup and deallocate resources
    subroutine cleanup_and_exit(xc_func, V_ext, results)
        type(xc_lsda_t), intent(inout) :: xc_func
        real(dp), allocatable, intent(inout) :: V_ext(:)
        type(scf_results_t), intent(inout) :: results
        integer :: ierr
        
        call cleanup_scf_results(results, ierr)
        call xc_lsda_destroy(xc_func)
        if (allocated(V_ext)) deallocate(V_ext)
    end subroutine cleanup_and_exit

end program lsdaks
