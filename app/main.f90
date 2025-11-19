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
    call check_xc_table(sys_params%U, table_file, table_exists)
    
    if (.not. table_exists) then
        print '(A)', ""
        print '(A)', "=========================================="
        print '(A)', "ERROR: XC table not found!"
        print '(A)', "=========================================="
        print '(A,F0.2)', "  Requested U = ", sys_params%U
        print '(A)', ""
        print '(A)', "Please generate the XC table first using:"
        print '(A,F0.2)', "  fpm run generate_xc_table -- --U ", sys_params%U
        print '(A)', ""
        print '(A)', "Available tables are stored in:"
        print '(A)', "  data/tables/fortran_native/"
        print '(A)', ""
        call print_available_tables()
        stop 1
    end if
    
    print '(A)', ""
    print '(A,A)', "Loading XC table: ", trim(table_file)
    
    call xc_lsda_init(xc_func, table_file, ierr)
    if (ierr /= ERROR_SUCCESS) then
        print *, "ERROR: Failed to initialize XC functional"
        stop 1
    end if
    
    print '(A)', "  ✓ XC functional initialized"

    allocate(V_ext(sys_params%L))

    ! Prepare potential parameters based on type
    select case (trim(inputs%potential_type))
    case ("uniform")
        allocate(pot_params(1))
        pot_params(1) = inputs%V0
    case ("harmonic")
        allocate(pot_params(1))
        pot_params(1) = inputs%V0  ! Using V0 as spring constant k
    case default
        ! For other potentials, allocate minimal params
        allocate(pot_params(1))
        pot_params(1) = inputs%V0
    end select

    seed = -1  ! Use system time for random potentials

    call create_potential(inputs%potential_type, pot_params, sys_params%L, seed, V_ext, ierr)

    deallocate(pot_params)

    if (ierr /= ERROR_SUCCESS) then
        print *, "ERROR: Failed to create external potential"
        call xc_lsda_destroy(xc_func)
        deallocate(V_ext)
        stop 1
    end if
    
    print '(A,A)', "  ✓ External potential created: ", trim(inputs%potential_type)
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
        print '(A)', ""
        
        print '(A)', "SCF Parameters:"
        print '(A,I0)', "  Max iterations:   ", scf_params%max_iter
        print '(A,ES10.2)', "  Density tol:      ", scf_params%density_tol
        print '(A,ES10.2)', "  Energy tol:       ", scf_params%energy_tol
        print '(A,F0.3)', "  Mixing alpha:     ", scf_params%mixing_alpha
        print '(A,L1)', "  Verbose:          ", scf_params%verbose
        print '(A)', "----------------------------------------"
    end subroutine print_configuration

    !> Check if XC table exists for given U
    subroutine check_xc_table(U, table_file, exists)
        real(dp), intent(in) :: U
        character(len=*), intent(out) :: table_file
        logical, intent(out) :: exists
        
        write(table_file, '(A,F0.2,A)') 'data/tables/fortran_native/xc_table_u', U, '.dat'
        
        inquire(file=table_file, exist=exists)
    end subroutine check_xc_table

    !> Print list of available XC tables
    subroutine print_available_tables()
        real(dp), dimension(25) :: available_U
        integer :: i
        logical :: exists
        character(len=256) :: filename
        
        available_U = [1.00_dp, 1.10_dp, 2.00_dp, 3.00_dp, 4.00_dp, &
                      4.10_dp, 5.00_dp, 5.90_dp, 6.00_dp, 6.10_dp, &
                      6.90_dp, 7.00_dp, 7.10_dp, 7.90_dp, 8.00_dp, &
                      8.10_dp, 8.90_dp, 9.00_dp, 9.10_dp, 10.00_dp, &
                      12.00_dp, 14.00_dp, 16.00_dp, 18.00_dp, 20.00_dp]
        
        print '(A)', "Common U values with available tables:"
        print '(A)', ""
        
        do i = 1, 25
            write(filename, '(A,F0.2,A)') 'data/tables/fortran_native/xc_table_u', &
                                         available_U(i), '.dat'
            inquire(file=filename, exist=exists)
            
            if (exists) then
                if (mod(i-1, 5) == 0) print '(A)', ""
                write(*, '(A,F5.2)', advance='no') "  U = ", available_U(i)
            end if
        end do
        
        print '(A)', ""
        print '(A)', ""
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