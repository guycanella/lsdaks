!> Module for writing simulation results to files
!!
!! Provides formatted output for:
!! - Density profiles (spin-up, spin-down, total)
!! - Eigenvalues and eigenvectors
!! - Convergence history
!! - Summary information
module output_writer
    use lsda_constants, only: dp
    use lsda_types, only: system_params_t
    use kohn_sham_cycle, only: scf_results_t, scf_params_t
    use input_parser, only: input_params_t
    use lsda_errors, only: ERROR_SUCCESS, ERROR_FILE_WRITE
    implicit none
    private

    public :: write_results
    public :: write_density_profile
    public :: write_eigenvalues
    public :: write_convergence_history
    public :: write_summary

contains

    !> Write all simulation results
    !!
    !! @param[in] results      SCF results
    !! @param[in] sys_params   System parameters
    !! @param[in] inputs       Input parameters
    !! @param[out] ierr        Error code
    subroutine write_results(results, sys_params, inputs, ierr)
        type(scf_results_t), intent(in) :: results
        type(system_params_t), intent(in) :: sys_params
        type(input_params_t), intent(in) :: inputs
        integer, intent(out) :: ierr
        
        character(len=256) :: prefix
        
        ierr = ERROR_SUCCESS
        prefix = trim(inputs%output_prefix)
        
        call write_summary(results, sys_params, inputs, ierr)
        if (ierr /= ERROR_SUCCESS) return
        
        if (inputs%save_density) then
            call write_density_profile(results, sys_params, prefix, ierr)
            if (ierr /= ERROR_SUCCESS) return
        end if
        
        if (inputs%save_eigenvalues) then
            call write_eigenvalues(results, sys_params, prefix, ierr)
            if (ierr /= ERROR_SUCCESS) return
        end if
        
        if (inputs%store_history .and. allocated(results%history%density_norms)) then
            call write_convergence_history(results, prefix, ierr)
            if (ierr /= ERROR_SUCCESS) return
        end if
        
    end subroutine write_results

    !> Write summary to stdout and file
    !!
    !! @param[in] results      SCF results
    !! @param[in] sys_params   System parameters
    !! @param[in] inputs       Input parameters
    !! @param[out] ierr        Error code
    subroutine write_summary(results, sys_params, inputs, ierr)
        type(scf_results_t), intent(in) :: results
        type(system_params_t), intent(in) :: sys_params
        type(input_params_t), intent(in) :: inputs
        integer, intent(out) :: ierr

        character(len=256) :: filename
        integer :: io_unit, io_stat
        real(dp) :: total_density
        
        ierr = ERROR_SUCCESS
        
        print '(A)', ""
        print '(A)', "=========================================="
        print '(A)', "         SIMULATION RESULTS"
        print '(A)', "=========================================="
        print '(A)', ""
        
        print '(A)', "System Parameters:"
        print '(A,I0)', "  L (sites):        ", sys_params%L
        print '(A,I0)', "  N_up:             ", sys_params%Nup
        print '(A,I0)', "  N_down:           ", sys_params%Ndown
        print '(A,I0)', "  N_total:          ", sys_params%Nup + sys_params%Ndown
        print '(A,F0.4)', "  U:                ", sys_params%U
        print '(A,I0)', "  BC:               ", sys_params%bc
        if (sys_params%bc == 2) then
            print '(A,F0.4)', "  Phase:            ", sys_params%phase
        end if
        print '(A)', ""
        
        print '(A)', "SCF Convergence:"
        if (results%converged) then
            print '(A)', "  Status:           ✓ CONVERGED"
        else
            print '(A)', "  Status:           ✗ NOT CONVERGED"
        end if
        print '(A,I0)', "  Iterations:       ", results%n_iterations
        print '(A,ES12.4)', "  Final |Δn|:       ", results%final_density_error
        print '(A,F16.8)', "  Final Energy:     ", results%final_energy
        print '(A)', ""
        
        if (allocated(results%density_up)) then
            total_density = sum(results%density_up) + sum(results%density_down)
            print '(A)', "Density Check:"
            print '(A,F12.6)', "  ∫n_up dx:         ", sum(results%density_up)
            print '(A,F12.6)', "  ∫n_down dx:       ", sum(results%density_down)
            print '(A,F12.6)', "  ∫n_total dx:      ", total_density
            print '(A,F12.6)', "  Expected N:       ", real(sys_params%Nup + sys_params%Ndown, dp)
            print '(A,ES12.4)', "  Error:            ", abs(total_density - real(sys_params%Nup + sys_params%Ndown, dp))
            print '(A)', ""
        end if
        
        print '(A)', "Output Files:"
        print '(A,A)', "  Prefix:           ", trim(inputs%output_prefix)
        if (inputs%save_density) then
            print '(A,A)', "  Density:          ", trim(inputs%output_prefix) // "_density.dat"
        end if
        if (inputs%save_eigenvalues) then
            print '(A,A)', "  Eigenvalues:      ", trim(inputs%output_prefix) // "_eigenvalues.dat"
        end if
        if (inputs%store_history) then
            print '(A,A)', "  Convergence:      ", trim(inputs%output_prefix) // "_convergence.dat"
        end if
        print '(A)', ""
        print '(A)', "=========================================="
        print '(A)', ""
        
        filename = trim(inputs%output_prefix) // "_summary.txt"
        open(newunit=io_unit, file=filename, status='replace', iostat=io_stat)
        if (io_stat /= 0) then
            print *, "WARNING: Could not write summary file: ", trim(filename)
            ierr = ERROR_FILE_WRITE
            return
        end if
        
        write(io_unit, '(A)') "LSDA-Hubbard Simulation Summary"
        write(io_unit, '(A)') "================================"
        write(io_unit, '(A)') ""
        write(io_unit, '(A,I0)') "L = ", sys_params%L
        write(io_unit, '(A,I0)') "Nup = ", sys_params%Nup
        write(io_unit, '(A,I0)') "Ndown = ", sys_params%Ndown
        write(io_unit, '(A,F0.4)') "U = ", sys_params%U
        write(io_unit, '(A)') ""
        
        if (results%converged) then
            write(io_unit, '(A)') "SCF: CONVERGED"
        else
            write(io_unit, '(A)') "SCF: NOT CONVERGED"
        end if
        write(io_unit, '(A,I0)') "Iterations: ", results%n_iterations
        write(io_unit, '(A,ES12.4)') "Final |Δn|: ", results%final_density_error
        write(io_unit, '(A,F16.8)') "Final Energy: ", results%final_energy
        
        close(io_unit)
        
        ierr = ERROR_SUCCESS
        
    end subroutine write_summary

    !> Write density profile to file
    !!
    !! Format: site n_up(i) n_down(i) n_total(i)
    !!
    !! @param[in] results      SCF results
    !! @param[in] sys_params   System parameters
    !! @param[in] prefix       Output file prefix
    !! @param[out] ierr        Error code
    subroutine write_density_profile(results, sys_params, prefix, ierr)
        type(scf_results_t), intent(in) :: results
        type(system_params_t), intent(in) :: sys_params
        character(len=*), intent(in) :: prefix
        integer, intent(out) :: ierr

        character(len=256) :: filename
        integer :: io_unit, io_stat, i
        
        ierr = ERROR_SUCCESS
        
        if (.not. allocated(results%density_up)) then
            print *, "WARNING: Density not available, skipping density output"
            return
        end if
        
        filename = trim(prefix) // "_density.dat"

        open(newunit=io_unit, file=filename, status='replace', iostat=io_stat)
        if (io_stat /= 0) then
            print *, "ERROR: Could not write density file: ", trim(filename)
            ierr = ERROR_FILE_WRITE
            return
        end if
        
        write(io_unit, '(A)') "# Density profile from LSDA-Hubbard calculation"
        write(io_unit, '(A,I0)') "# L = ", sys_params%L
        write(io_unit, '(A,I0)') "# Nup = ", sys_params%Nup
        write(io_unit, '(A,I0)') "# Ndown = ", sys_params%Ndown
        write(io_unit, '(A,F0.4)') "# U = ", sys_params%U
        write(io_unit, '(A)') "#"
        write(io_unit, '(A)') "# Columns: site  n_up  n_down  n_total"
        
        do i = 1, sys_params%L
            write(io_unit, '(I6,3ES20.10)') i, results%density_up(i), &
                                             results%density_down(i), &
                                             results%density_up(i) + results%density_down(i)
        end do
        
        close(io_unit)
        
        print '(A,A)', "  Density profile written to: ", trim(filename)
        
        ierr = ERROR_SUCCESS
        
    end subroutine write_density_profile

    !> Write eigenvalues to file
    !!
    !! Format: index spin eigenvalue
    !!
    !! @param[in] results      SCF results
    !! @param[in] sys_params   System parameters
    !! @param[in] prefix       Output file prefix
    !! @param[out] ierr        Error code
    subroutine write_eigenvalues(results, sys_params, prefix, ierr)
        type(scf_results_t), intent(in) :: results
        type(system_params_t), intent(in) :: sys_params
        character(len=*), intent(in) :: prefix
        integer, intent(out) :: ierr

        character(len=256) :: filename
        integer :: io_unit, io_stat, i, L
        
        ierr = ERROR_SUCCESS
        
        if (.not. allocated(results%eigvals)) then
            print *, "WARNING: Eigenvalues not available, skipping eigenvalue output"
            return
        end if
        
        L = sys_params%L
        filename = trim(prefix) // "_eigenvalues.dat"

        open(newunit=io_unit, file=filename, status='replace', iostat=io_stat)
        if (io_stat /= 0) then
            print *, "ERROR: Could not write eigenvalues file: ", trim(filename)
            ierr = ERROR_FILE_WRITE
            return
        end if
        
        ! Write header
        write(io_unit, '(A)') "# Eigenvalues from LSDA-Hubbard calculation"
        write(io_unit, '(A,I0)') "# L = ", L
        write(io_unit, '(A,I0)') "# Nup = ", sys_params%Nup
        write(io_unit, '(A,I0)') "# Ndown = ", sys_params%Ndown
        write(io_unit, '(A,F0.4)') "# U = ", sys_params%U
        write(io_unit, '(A)') "#"
        write(io_unit, '(A)') "# First L eigenvalues: spin-up"
        write(io_unit, '(A)') "# Last L eigenvalues: spin-down"
        write(io_unit, '(A)') "#"
        write(io_unit, '(A)') "# Columns: index  spin  eigenvalue  occupied"
        
        do i = 1, L
            if (i <= sys_params%Nup) then
                write(io_unit, '(I6,A8,ES20.10,A8)') i, "up", results%eigvals(i), "yes"
            else
                write(io_unit, '(I6,A8,ES20.10,A8)') i, "up", results%eigvals(i), "no"
            end if
        end do
        
        do i = 1, L
            if (i <= sys_params%Ndown) then
                write(io_unit, '(I6,A8,ES20.10,A8)') i, "down", results%eigvals(L+i), "yes"
            else
                write(io_unit, '(I6,A8,ES20.10,A8)') i, "down", results%eigvals(L+i), "no"
            end if
        end do
        
        close(io_unit)
        
        print '(A,A)', "  Eigenvalues written to: ", trim(filename)
        
        ierr = ERROR_SUCCESS
        
    end subroutine write_eigenvalues

    !> Write convergence history to file
    !!
    !! Format: iteration density_error energy
    !!
    !! @param[in] results  SCF results
    !! @param[in] prefix   Output file prefix
    !! @param[out] ierr    Error code
    subroutine write_convergence_history(results, prefix, ierr)
        type(scf_results_t), intent(in) :: results
        character(len=*), intent(in) :: prefix
        integer, intent(out) :: ierr

        character(len=256) :: filename
        integer :: io_unit, io_stat, i
        
        ierr = ERROR_SUCCESS
        
        if (.not. allocated(results%history%density_norms)) then
            return
        end if
        
        filename = trim(prefix) // "_convergence.dat"

        open(newunit=io_unit, file=filename, status='replace', iostat=io_stat)
        if (io_stat /= 0) then
            print *, "ERROR: Could not write convergence file: ", trim(filename)
            ierr = ERROR_FILE_WRITE
            return
        end if
        
        write(io_unit, '(A)') "# SCF convergence history"
        write(io_unit, '(A)') "#"
        write(io_unit, '(A)') "# Columns: iteration  |Δn|  energy"
        
        do i = 1, results%history%current_iter
            write(io_unit, '(I6,2ES20.10)') i, &
                                            results%history%density_norms(i), &
                                            results%history%energies(i)
        end do
        
        close(io_unit)
        
        print '(A,A)', "  Convergence history written to: ", trim(filename)
        
        ierr = ERROR_SUCCESS
        
    end subroutine write_convergence_history

end module output_writer