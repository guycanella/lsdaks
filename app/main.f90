!> Main program for LSDA-Hubbard calculations
!!
!! LSDAKS: Local Spin Density Approximation - Kohn-Sham solver
!! for the 1D Hubbard model using Bethe Ansatz-based XC functionals.
program lsdaks
    use lsda_constants, only: dp, PI
    use lsda_types, only: system_params_t
    use lsda_errors, only: ERROR_SUCCESS, ERROR_FILE_NOT_FOUND, ERROR_CONVERGENCE_FAILED
    use input_parser
    use output_writer
    use kohn_sham_cycle
    use xc_lsda, only: xc_lsda_t, xc_lsda_init, xc_lsda_destroy
    use potential_factory, only: create_potential
    use potential_impurity, only: potential_impurity_single, potential_impurity_random
    use boundary_conditions, only: BC_OPEN, BC_PERIODIC, BC_TWISTED
    implicit none
    
    ! Main variables
    type(input_params_t) :: inputs
    type(system_params_t) :: sys_params
    type(scf_params_t) :: scf_params
    type(scf_results_t) :: results
    type(xc_lsda_t) :: xc_func
    
    real(dp), allocatable :: V_ext(:)
    real(dp), allocatable :: pot_params(:)
    integer, allocatable :: imp_positions(:)
    character(len=256) :: table_file
    integer :: ierr, seed
    logical :: table_exists
    
    call print_banner()
    
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
    
    call convert_to_system_params(inputs, sys_params, ierr)
    if (ierr /= ERROR_SUCCESS) then
        print *, "ERROR: Failed to convert system parameters"
        stop 1
    end if
    
    call convert_to_scf_params(inputs, scf_params)
    call print_configuration(inputs, sys_params, scf_params)
    
    ! Check XC table using |U| (tables are symmetric)
    call check_xc_table(abs(sys_params%U), table_file, table_exists)
    
    if (.not. table_exists) then
        print '(A)', ""
        print '(A)', "=========================================="
        print '(A)', "ERROR: XC table not found!"
        print '(A)', "=========================================="
        print '(A,F0.2)', "  Requested |U| = ", abs(sys_params%U)
        print '(A)', ""
        print '(A)', "Please generate the XC table first using:"
        print '(A,F0.2)', "  fpm run generate_xc_table -- --U ", abs(sys_params%U)
        print '(A)', ""
        call print_available_tables()
        stop 1
    end if
    
    print '(A)', ""
    print '(A,A)', "Loading XC table: ", trim(table_file)
    if (sys_params%U < 0.0_dp) then
        print '(A)', "  Note: Using table for |U| (attractive interaction)"
    end if
    
    call xc_lsda_init(xc_func, table_file, ierr)
    if (ierr /= ERROR_SUCCESS) then
        print *, "ERROR: Failed to initialize XC functional"
        stop 1
    end if
    
    print '(A)', "  ✓ XC functional initialized"

    allocate(V_ext(sys_params%L))

    ! Create external potential - special handling for impurities
    select case (trim(inputs%potential_type))
    
    case ('impurity')
        ! Random impurities with concentration
        seed = inputs%pot_seed
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
    
    case default
        ! Use factory for other potential types
        select case (trim(inputs%potential_type))
        case ("uniform")
            allocate(pot_params(1))
            pot_params(1) = inputs%V0
        case ("harmonic")
            allocate(pot_params(1))
            pot_params(1) = inputs%V0
        case ("random_uniform", "random_gaussian")
            allocate(pot_params(1))
            pot_params(1) = inputs%V0
        case ("barrier_single")
            allocate(pot_params(3))
            pot_params(1) = inputs%V0
            pot_params(2) = inputs%pot_center - inputs%pot_width/2.0_dp
            pot_params(3) = inputs%pot_center + inputs%pot_width/2.0_dp
        case default
            allocate(pot_params(1))
            pot_params(1) = inputs%V0
        end select

        seed = inputs%pot_seed

        call create_potential(inputs%potential_type, pot_params, sys_params%L, seed, V_ext, ierr)
        deallocate(pot_params)

        if (ierr /= ERROR_SUCCESS) then
            print *, "ERROR: Failed to create external potential"
            call xc_lsda_destroy(xc_func)
            deallocate(V_ext)
            stop 1
        end if
        
        print '(A,A)', "  ✓ External potential created: ", trim(inputs%potential_type)
    end select
    
    print '(A)', ""
    
    call init_scf_results(results, sys_params%L, scf_params%store_history, &
                         scf_params%max_iter, ierr)
    if (ierr /= ERROR_SUCCESS) then
        print *, "ERROR: Failed to initialize SCF results"
        call xc_lsda_destroy(xc_func)
        deallocate(V_ext)
        stop 1
    end if
    
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
    
    if (ierr /= ERROR_SUCCESS) then
        if (ierr == ERROR_CONVERGENCE_FAILED) then
            print '(A)', "=========================================="
            print '(A)', "WARNING: SCF did not converge!"
            print '(A)', "=========================================="
            print '(A,I0)', "  Iterations performed: ", results%n_iterations
            print '(A,ES12.4)', "  Final density error:  ", results%final_density_error
            print '(A,F16.8)', "  Final energy:         ", results%final_energy
            print '(A)', ""
            print '(A)', "Results may be unreliable."
            print '(A)', "Consider:"
            print '(A)', "  - Increasing max_iter"
            print '(A)', "  - Adjusting mixing_alpha"
            print '(A)', "  - Checking system parameters"
            print '(A)', ""
        else
            print *, "ERROR: SCF calculation failed with error code:", ierr
            call cleanup_and_exit(xc_func, V_ext, results)
            stop 1
        end if
    end if
    
    print '(A)', "=========================================="
    print '(A)', "  Writing Results"
    print '(A)', "=========================================="
    print '(A)', ""
    
    call write_results(results, sys_params, inputs, ierr)
    if (ierr /= ERROR_SUCCESS) then
        print *, "WARNING: Some output files could not be written"
    end if
    
    call cleanup_and_exit(xc_func, V_ext, results)
    
    print '(A)', ""
    print '(A)', "=========================================="
    print '(A)', "  Calculation Complete!"
    print '(A)', "=========================================="
    print '(A)', ""
    
contains

    subroutine print_banner()
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
        print '(A,F0.4)', "  Filling:          ", real(sys_params%Nup + sys_params%Ndown, dp) / real(sys_params%L, dp)
        print '(A,F0.4)', "  U:                ", sys_params%U
        if (sys_params%U < 0.0_dp) then
            print '(A)', "    (attractive interaction - pairing)"
        end if
        print '(A,A)', "  BC:               ", trim(inputs%bc_type)
        if (sys_params%bc == BC_TWISTED) then
            print '(A,F0.4,A)', "  Phase:            ", sys_params%phase, " π"
        end if
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
        print '(A,ES10.2)', "  Density tol:      ", scf_params%density_tol
        print '(A,ES10.2)', "  Energy tol:       ", scf_params%energy_tol
        print '(A,F0.3)', "  Mixing alpha:     ", scf_params%mixing_alpha
        print '(A,L1)', "  Verbose:          ", scf_params%verbose
        print '(A)', "----------------------------------------"
    end subroutine print_configuration

    !> Check if XC table exists for given U (uses |U|)
    subroutine check_xc_table(U, table_file, exists)
        real(dp), intent(in) :: U
        character(len=*), intent(out) :: table_file
        logical, intent(out) :: exists
        
        write(table_file, '(A,F0.2,A)') 'data/tables/fortran_native/xc_table_u', abs(U), '.dat'
        
        inquire(file=table_file, exist=exists)
    end subroutine check_xc_table

    !> Print list of ACTUALLY available XC tables (scan directory)
    subroutine print_available_tables()
        character(len=256) :: table_dir, filename, command
        integer :: io_stat, io_unit, i, n_found
        real(dp) :: U_value
        real(dp), allocatable :: found_U(:)
        logical :: exists
        character(len=10) :: U_str
        
        table_dir = 'data/tables/fortran_native/'
        
        inquire(file=trim(table_dir), exist=exists)
        if (.not. exists) then
            print '(A)', "Table directory not found: ", trim(table_dir)
            print '(A)', ""
            print '(A)', "Please create the directory and generate tables using:"
            print '(A)', "  mkdir -p data/tables/fortran_native"
            print '(A)', "  fpm run generate_xc_table -- --U <value>"
            return
        end if
        
        print '(A)', "Scanning for available XC tables in:"
        print '(A,A)', "  ", trim(table_dir)
        print '(A)', ""
        
        allocate(found_U(100))  ! Max 100 tables
        n_found = 0
        
        do i = 1, 20
            U_value = real(i, dp)
            write(filename, '(A,A,F0.2,A)') trim(table_dir), 'xc_table_u', U_value, '.dat'
            inquire(file=filename, exist=exists)
            if (exists) then
                n_found = n_found + 1
                found_U(n_found) = U_value
            end if
        end do
        
        do i = 1, 20
            U_value = real(i, dp) + 0.1_dp
            write(filename, '(A,A,F0.2,A)') trim(table_dir), 'xc_table_u', U_value, '.dat'
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