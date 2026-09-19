program test_output_writer
    use fortuno_serial, only: execute_serial_cmd_app
    implicit none

    call execute_serial_cmd_app(get_output_writer_tests())

contains

    function get_output_writer_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests

        tests = test_list([ &
            test("write_summary_success", test_write_summary_success), &
            test("write_summary_converged", test_write_summary_converged), &
            test("write_summary_not_converged", test_write_summary_not_converged), &
            test("write_density_profile", test_write_density_profile), &
            test("write_density_not_allocated", test_write_density_not_allocated), &
            test("write_eigenvalues", test_write_eigenvalues), &
            test("write_eigenvalues_only_computed_levels", test_write_eigenvalues_only_computed_levels), &
            test("write_eigenvalues_not_allocated", test_write_eigenvalues_not_allocated), &
            test("write_convergence_history", test_write_convergence_history), &
            test("write_convergence_not_allocated", test_write_convergence_not_allocated), &
            test("write_results_all_enabled", test_write_results_all_enabled), &
            test("write_results_density_disabled", test_write_results_density_disabled), &
            test("density_profile_format", test_density_profile_format), &
            test("eigenvalues_format_occupied", test_eigenvalues_format_occupied), &
            test("eigenvalues_format_unoccupied", test_eigenvalues_format_unoccupied), &
            test("outputs_record_smoothing_and_tolerance", &
                 test_outputs_record_smoothing_and_tolerance), &
            test("outputs_record_external_potential_provenance", &
                 test_outputs_record_external_potential_provenance), &
            test("provenance_records_type_specific_params", &
                 test_provenance_records_type_specific_params), &
            test("quasiperiodic_provenance_records_generator_params", &
                 test_quasiperiodic_provenance_records_generator_params), &
            test("convergence_header_identifies_system", &
                 test_convergence_header_identifies_system), &
            test("random_uniform_provenance_states_no_false_distribution", &
                 test_random_uniform_no_false_distribution), &
            test("not_converged_run_writes_all_files", &
                 test_not_converged_run_writes_all_files), &
            test("energy_per_site_uses_system_size", &
                 test_energy_per_site_uses_system_size) &
        ])
    end function get_output_writer_tests


    !> Does `filename` contain a line holding `needle`?
    !!
    !! @param[in] filename Path of the file to scan
    !! @param[in] needle   Substring to look for
    !! @return             .true. if some line of the file contains `needle`
    function file_contains(filename, needle) result(found)
        character(len=*), intent(in) :: filename, needle
        logical :: found

        integer :: io_unit, io_stat
        character(len=512) :: line

        found = .false.
        open(newunit=io_unit, file=filename, status='old', iostat=io_stat)
        if (io_stat /= 0) return

        do
            read(io_unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) exit
            if (index(line, needle) > 0) found = .true.
        end do

        close(io_unit)
    end function file_contains


    !> Delete a file if it exists (test housekeeping)
    !!
    !! @param[in] filename Path of the file to remove
    subroutine remove_file(filename)
        character(len=*), intent(in) :: filename
        integer :: io_unit, io_stat

        open(newunit=io_unit, file=filename, status='old', iostat=io_stat)
        if (io_stat == 0) close(io_unit, status='delete')
    end subroutine remove_file


    !> Test write_summary creates file successfully
    subroutine test_write_summary_success()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use lsda_types, only: system_params_t
        use input_parser, only: input_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        type(input_params_t) :: inputs
        integer :: ierr
        logical :: file_exists

        ! Setup minimal data
        sys_params%L = 10
        sys_params%Nup = 5
        sys_params%Ndown = 5
        sys_params%U = 4.0_dp
        sys_params%bc = 1

        results%converged = .true.
        results%n_iterations = 25
        results%final_density_error = 1.0e-9_dp
        results%final_energy = -12.5_dp

        inputs%output_prefix = 'test_summary'
        inputs%save_density = .false.
        inputs%save_eigenvalues = .false.
        inputs%store_history = .false.

        call write_summary(results, sys_params, inputs, ierr)

        call check(ierr == ERROR_SUCCESS, "write_summary should succeed")

        inquire(file='test_summary_summary.txt', exist=file_exists)
        call check(file_exists, "Summary file should exist")

        ! Cleanup
        if (file_exists) then
            open(unit=99, file='test_summary_summary.txt', status='old')
            close(99, status='delete')
        end if
    end subroutine test_write_summary_success


    !> Test write_summary with converged results
    subroutine test_write_summary_converged()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use lsda_types, only: system_params_t
        use input_parser, only: input_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        type(input_params_t) :: inputs
        integer :: ierr
        logical :: file_exists

        sys_params%L = 5
        sys_params%Nup = 3
        sys_params%Ndown = 2
        sys_params%U = 2.0_dp
        sys_params%bc = 1

        allocate(results%density_up(5))
        allocate(results%density_down(5))

        results%density_up = [0.6_dp, 0.6_dp, 0.6_dp, 0.6_dp, 0.6_dp]
        results%density_down = [0.4_dp, 0.4_dp, 0.4_dp, 0.4_dp, 0.4_dp]
        results%converged = .true.
        results%n_iterations = 15
        results%final_density_error = 1.0e-10_dp
        results%final_energy = -8.0_dp

        inputs%output_prefix = 'test_conv'
        inputs%save_density = .false.
        inputs%save_eigenvalues = .false.
        inputs%store_history = .false.

        call write_summary(results, sys_params, inputs, ierr)

        call check(ierr == ERROR_SUCCESS, "Converged summary should succeed")

        inquire(file='test_conv_summary.txt', exist=file_exists)
        call check(file_exists, "Converged summary file should exist")

        ! Cleanup
        if (file_exists) then
            open(unit=99, file='test_conv_summary.txt', status='old')
            close(99, status='delete')
        end if

        deallocate(results%density_up)
        deallocate(results%density_down)
    end subroutine test_write_summary_converged


    !> Test write_summary with non-converged results
    subroutine test_write_summary_not_converged()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use lsda_types, only: system_params_t
        use input_parser, only: input_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        type(input_params_t) :: inputs
        integer :: ierr
        logical :: file_exists

        sys_params%L = 8
        sys_params%Nup = 4
        sys_params%Ndown = 4
        sys_params%U = 6.0_dp
        sys_params%bc = 1

        results%converged = .false.  ! Not converged
        results%n_iterations = 100
        results%final_density_error = 1.0e-5_dp
        results%final_energy = -10.0_dp

        inputs%output_prefix = 'test_noconv'
        inputs%save_density = .false.
        inputs%save_eigenvalues = .false.
        inputs%store_history = .false.

        call write_summary(results, sys_params, inputs, ierr)

        call check(ierr == ERROR_SUCCESS, "Non-converged summary should succeed")

        inquire(file='test_noconv_summary.txt', exist=file_exists)
        call check(file_exists, "Non-converged summary file should exist")

        ! Cleanup
        if (file_exists) then
            open(unit=99, file='test_noconv_summary.txt', status='old')
            close(99, status='delete')
        end if
    end subroutine test_write_summary_not_converged


    !> Test write_density_profile creates file with correct format
    subroutine test_write_density_profile()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use lsda_types, only: system_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        integer :: ierr
        logical :: file_exists

        sys_params%L = 5
        sys_params%Nup = 3
        sys_params%Ndown = 2
        sys_params%U = 4.0_dp

        allocate(results%density_up(5))
        allocate(results%density_down(5))

        results%density_up = [0.6_dp, 0.6_dp, 0.6_dp, 0.6_dp, 0.6_dp]
        results%density_down = [0.4_dp, 0.4_dp, 0.4_dp, 0.4_dp, 0.4_dp]

        call write_density_profile(results, sys_params, 'test_dens', ierr)

        call check(ierr == ERROR_SUCCESS, "write_density_profile should succeed")

        inquire(file='test_dens_density.dat', exist=file_exists)
        call check(file_exists, "Density file should exist")

        ! Cleanup
        if (file_exists) then
            open(unit=99, file='test_dens_density.dat', status='old')
            close(99, status='delete')
        end if

        deallocate(results%density_up)
        deallocate(results%density_down)
    end subroutine test_write_density_profile


    !> Test write_density_profile when density not allocated
    subroutine test_write_density_not_allocated()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use lsda_types, only: system_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        integer :: ierr

        sys_params%L = 5
        sys_params%Nup = 3
        sys_params%Ndown = 2
        sys_params%U = 4.0_dp

        ! Don't allocate density arrays

        call write_density_profile(results, sys_params, 'test_nodens', ierr)

        ! Should return successfully but without writing
        call check(ierr == ERROR_SUCCESS, "Should handle unallocated density gracefully")
    end subroutine test_write_density_not_allocated


    !> Test write_eigenvalues creates file
    subroutine test_write_eigenvalues()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use lsda_types, only: system_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        integer :: ierr
        logical :: file_exists

        sys_params%L = 4
        sys_params%Nup = 2
        sys_params%Ndown = 2
        sys_params%U = 4.0_dp

        allocate(results%eigvals(8))  ! 2*L
        results%eigvals = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, &
                          -1.5_dp, -0.5_dp, 0.5_dp, 1.5_dp]

        call write_eigenvalues(results, sys_params, 'test_eig', ierr)

        call check(ierr == ERROR_SUCCESS, "write_eigenvalues should succeed")

        inquire(file='test_eig_eigenvalues.dat', exist=file_exists)
        call check(file_exists, "Eigenvalues file should exist")

        ! Cleanup
        if (file_exists) then
            open(unit=99, file='test_eig_eigenvalues.dat', status='old')
            close(99, status='delete')
        end if

        deallocate(results%eigvals)
    end subroutine test_write_eigenvalues


    !> Regression: only eigenvalues actually computed by the SCF are written.
    !!
    !! A partial diagonalization must not make absent high-energy levels appear
    !! as zero-valued eigenvalues in the output file.
    subroutine test_write_eigenvalues_only_computed_levels()
        use fortuno_serial, only: check => serial_check
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        use output_writer, only: write_eigenvalues
        use lsda_types, only: system_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        integer :: ierr, io_unit, io_stat, n_data, i
        character(len=256) :: line

        sys_params%L = 20
        sys_params%Nup = 4
        sys_params%Ndown = 4
        allocate(results%eigvals(40))
        results%eigvals = ieee_value(0.0_dp, ieee_quiet_nan)
        results%eigvals(1:5) = [(real(i, dp), i = 1, 5)]
        results%eigvals(21:25) = [(real(i, dp), i = 6, 10)]

        call write_eigenvalues(results, sys_params, 'test_eig_partial', ierr)
        call check(ierr == ERROR_SUCCESS, "partial eigenvalue output should succeed")

        n_data = 0
        open(newunit=io_unit, file='test_eig_partial_eigenvalues.dat', status='old', iostat=io_stat)
        call check(io_stat == 0, "partial eigenvalue file should open")
        if (io_stat == 0) then
            do
                read(io_unit, '(A)', iostat=io_stat) line
                if (io_stat /= 0) exit
                if (line(1:1) /= '#') n_data = n_data + 1
            end do
            close(io_unit, status='delete')
        end if
        call check(n_data == 10, "writer must emit exactly the computed eigenvalues")
        deallocate(results%eigvals)
    end subroutine test_write_eigenvalues_only_computed_levels


    !> Test write_eigenvalues when eigenvalues not allocated
    subroutine test_write_eigenvalues_not_allocated()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use lsda_types, only: system_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        integer :: ierr

        sys_params%L = 4
        sys_params%Nup = 2
        sys_params%Ndown = 2
        sys_params%U = 4.0_dp

        ! Don't allocate eigvals

        call write_eigenvalues(results, sys_params, 'test_noeig', ierr)

        ! Should return successfully but without writing
        call check(ierr == ERROR_SUCCESS, "Should handle unallocated eigenvalues gracefully")
    end subroutine test_write_eigenvalues_not_allocated


    !> Test write_convergence_history
    subroutine test_write_convergence_history()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        integer :: ierr
        logical :: file_exists

        allocate(results%history%density_norms(3))
        allocate(results%history%energies(3))

        results%history%current_iter = 3
        results%history%density_norms = [1.0e-2_dp, 1.0e-4_dp, 1.0e-6_dp]
        results%history%energies = [-10.0_dp, -12.0_dp, -12.5_dp]

        call write_convergence_history(results, 'test_conv_hist', ierr)

        call check(ierr == ERROR_SUCCESS, "write_convergence_history should succeed")

        inquire(file='test_conv_hist_convergence.dat', exist=file_exists)
        call check(file_exists, "Convergence file should exist")

        ! Cleanup
        if (file_exists) then
            open(unit=99, file='test_conv_hist_convergence.dat', status='old')
            close(99, status='delete')
        end if

        deallocate(results%history%density_norms)
        deallocate(results%history%energies)
    end subroutine test_write_convergence_history


    !> Test write_convergence_history when history not allocated
    subroutine test_write_convergence_not_allocated()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        integer :: ierr

        ! Don't allocate history

        call write_convergence_history(results, 'test_nohist', ierr)

        ! Should return successfully but without writing
        call check(ierr == ERROR_SUCCESS, "Should handle unallocated history gracefully")
    end subroutine test_write_convergence_not_allocated


    !> Test write_results with all outputs enabled
    subroutine test_write_results_all_enabled()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use lsda_types, only: system_params_t
        use input_parser, only: input_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        type(input_params_t) :: inputs
        integer :: ierr
        logical :: file1, file2, file3, file4

        ! Setup
        sys_params%L = 3
        sys_params%Nup = 2
        sys_params%Ndown = 1
        sys_params%U = 4.0_dp
        sys_params%bc = 1

        results%converged = .true.
        results%n_iterations = 10
        results%final_density_error = 1.0e-10_dp
        results%final_energy = -5.0_dp

        allocate(results%density_up(3))
        allocate(results%density_down(3))
        allocate(results%eigvals(6))
        allocate(results%history%density_norms(2))
        allocate(results%history%energies(2))

        results%density_up = [0.7_dp, 0.7_dp, 0.6_dp]
        results%density_down = [0.3_dp, 0.3_dp, 0.4_dp]
        results%eigvals = [-1.0_dp, 0.0_dp, 1.0_dp, -0.5_dp, 0.5_dp, 1.5_dp]
        results%history%current_iter = 2
        results%history%density_norms = [1.0e-5_dp, 1.0e-10_dp]
        results%history%energies = [-4.8_dp, -5.0_dp]

        inputs%output_prefix = 'test_all'
        inputs%save_density = .true.
        inputs%save_eigenvalues = .true.
        inputs%store_history = .true.

        call write_results(results, sys_params, inputs, ierr)

        call check(ierr == ERROR_SUCCESS, "write_results should succeed")

        inquire(file='test_all_summary.txt', exist=file1)
        inquire(file='test_all_density.dat', exist=file2)
        inquire(file='test_all_eigenvalues.dat', exist=file3)
        inquire(file='test_all_convergence.dat', exist=file4)

        call check(file1, "Summary should exist")
        call check(file2, "Density should exist")
        call check(file3, "Eigenvalues should exist")
        call check(file4, "Convergence should exist")

        ! Cleanup
        if (file1) then
            open(unit=91, file='test_all_summary.txt', status='old')
            close(91, status='delete')
        end if
        if (file2) then
            open(unit=92, file='test_all_density.dat', status='old')
            close(92, status='delete')
        end if
        if (file3) then
            open(unit=93, file='test_all_eigenvalues.dat', status='old')
            close(93, status='delete')
        end if
        if (file4) then
            open(unit=94, file='test_all_convergence.dat', status='old')
            close(94, status='delete')
        end if

        deallocate(results%density_up)
        deallocate(results%density_down)
        deallocate(results%eigvals)
        deallocate(results%history%density_norms)
        deallocate(results%history%energies)
    end subroutine test_write_results_all_enabled


    !> Test write_results with density disabled
    subroutine test_write_results_density_disabled()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use lsda_types, only: system_params_t
        use input_parser, only: input_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        type(input_params_t) :: inputs
        integer :: ierr
        logical :: file_exists

        sys_params%L = 3
        sys_params%Nup = 2
        sys_params%Ndown = 1
        sys_params%U = 4.0_dp
        sys_params%bc = 1

        results%converged = .true.
        results%n_iterations = 10
        results%final_density_error = 1.0e-10_dp
        results%final_energy = -5.0_dp

        inputs%output_prefix = 'test_nodens_out'
        inputs%save_density = .false.  ! Disabled
        inputs%save_eigenvalues = .false.
        inputs%store_history = .false.

        call write_results(results, sys_params, inputs, ierr)

        call check(ierr == ERROR_SUCCESS, "write_results should succeed")

        ! Density file should NOT exist
        inquire(file='test_nodens_out_density.dat', exist=file_exists)
        call check(.not. file_exists, "Density file should not exist when disabled")

        ! Cleanup summary
        inquire(file='test_nodens_out_summary.txt', exist=file_exists)
        if (file_exists) then
            open(unit=99, file='test_nodens_out_summary.txt', status='old')
            close(99, status='delete')
        end if
    end subroutine test_write_results_density_disabled


    !> REGRESSION (T9): a non-convergent run gets the SAME files as a convergent one
    !!
    !! The failure path of the SCF cycle used to leave results%density_up,
    !! density_down and eigvals UNALLOCATED while still filling final_energy, so
    !! write_density_profile and write_eigenvalues took their "not available"
    !! branch and produced nothing. The run that most needs its density and its
    !! spectrum inspected - the one that did not converge - was the only one that
    !! did not get them on disk.
    !!
    !! The summary must also state the outcome in a machine-readable field:
    !! "SCF: NOT CONVERGED" contains the word CONVERGED, so grepping for the
    !! convergent case matches the divergent one too; `converged = F` cannot.
    subroutine test_not_converged_run_writes_all_files()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use lsda_types, only: system_params_t
        use input_parser, only: input_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        type(input_params_t) :: inputs
        integer :: ierr
        logical :: has_summary, has_density, has_eigvals, has_history

        sys_params%L = 4
        sys_params%Nup = 2
        sys_params%Ndown = 2
        sys_params%U = 4.0_dp
        sys_params%bc = 1

        ! Exactly the state the SCF cycle returns together with
        ! ERROR_CONVERGENCE_FAILED: converged = .false., the iteration budget
        ! exhausted, and the last available density / spectrum / energy.
        results%converged = .false.
        results%n_iterations = 200
        results%final_density_error = 3.0e-2_dp
        results%final_potential_residual = 5.0e-2_dp
        results%final_energy = -8.0_dp

        allocate(results%density_up(4))
        allocate(results%density_down(4))
        allocate(results%eigvals(8))
        allocate(results%history%density_norms(2))
        allocate(results%history%energies(2))

        results%density_up = [0.6_dp, 0.4_dp, 0.4_dp, 0.6_dp]
        results%density_down = [0.4_dp, 0.6_dp, 0.6_dp, 0.4_dp]
        results%eigvals = [-2.0_dp, -1.0_dp, 1.0_dp, 2.0_dp, &
                           -2.0_dp, -1.0_dp, 1.0_dp, 2.0_dp]
        results%history%current_iter = 2
        results%history%density_norms = [5.0e-2_dp, 3.0e-2_dp]
        results%history%energies = [-7.9_dp, -8.0_dp]

        inputs%output_prefix = 'test_t9_notconv'
        inputs%save_density = .true.
        inputs%save_eigenvalues = .true.
        inputs%store_history = .true.

        call write_results(results, sys_params, inputs, ierr)

        call check(ierr == ERROR_SUCCESS, &
                   "write_results must succeed for a non-convergent run")

        inquire(file='test_t9_notconv_summary.txt', exist=has_summary)
        inquire(file='test_t9_notconv_density.dat', exist=has_density)
        inquire(file='test_t9_notconv_eigenvalues.dat', exist=has_eigvals)
        inquire(file='test_t9_notconv_convergence.dat', exist=has_history)

        call check(has_summary, "A non-convergent run must still write its summary")
        call check(has_density, "A non-convergent run must still write its density profile")
        call check(has_eigvals, "A non-convergent run must still write its eigenvalues")
        call check(has_history, "A non-convergent run must still write its history")

        if (has_summary) then
            call check(file_contains('test_t9_notconv_summary.txt', 'converged = F'), &
                       "The summary must carry a machine-readable converged = F field")
            call check(file_contains('test_t9_notconv_summary.txt', 'SCF: NOT CONVERGED'), &
                       "The summary must state NOT CONVERGED in words as well")
        end if

        call remove_file('test_t9_notconv_summary.txt')
        call remove_file('test_t9_notconv_density.dat')
        call remove_file('test_t9_notconv_eigenvalues.dat')
        call remove_file('test_t9_notconv_convergence.dat')

        deallocate(results%density_up)
        deallocate(results%density_down)
        deallocate(results%eigvals)
        deallocate(results%history%density_norms)
        deallocate(results%history%energies)
    end subroutine test_not_converged_run_writes_all_files


    !> REGRESSION (T9): the energy per site is the total energy divided by L
    !!
    !! `results%final_energy` is the TOTAL energy. The summary used to divide it
    !! by `size(results%density_up)` and, when that array was absent, printed the
    !! TOTAL energy verbatim under the label "Final Energy per site" - wrong by a
    !! factor of L, with nothing in the file to reveal it. In the summary FILE the
    !! label did not mention "per site" at all, so the same number was reported
    !! under two different meanings in two different files of the same run.
    !!
    !! The test sets L = 4 and E_total = -8, so the two numbers (-8 and -2) cannot
    !! be confused, and deliberately leaves density_up UNALLOCATED, which is the
    !! configuration in which the divisor used to disappear entirely.
    subroutine test_energy_per_site_uses_system_size()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use lsda_types, only: system_params_t
        use input_parser, only: input_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        type(input_params_t) :: inputs
        integer :: ierr
        logical :: has_summary

        sys_params%L = 4
        sys_params%Nup = 2
        sys_params%Ndown = 2
        sys_params%U = 4.0_dp
        sys_params%bc = 1

        results%converged = .true.
        results%n_iterations = 12
        results%final_density_error = 1.0e-10_dp
        results%final_potential_residual = 1.0e-10_dp
        results%final_energy = -8.0_dp
        ! density_up / density_down deliberately NOT allocated.

        inputs%output_prefix = 'test_t9_energy'
        inputs%save_density = .false.
        inputs%save_eigenvalues = .false.
        inputs%store_history = .false.

        call write_summary(results, sys_params, inputs, ierr)

        call check(ierr == ERROR_SUCCESS, "write_summary should succeed")

        inquire(file='test_t9_energy_summary.txt', exist=has_summary)
        call check(has_summary, "Summary file should exist")

        if (has_summary) then
            call check(file_contains('test_t9_energy_summary.txt', 'Final Total Energy:'), &
                       "The summary must report the total energy under its own label")
            call check(file_contains('test_t9_energy_summary.txt', '-8.000000000000'), &
                       "The reported total energy must be E_total = -8")
            call check(file_contains('test_t9_energy_summary.txt', 'Final Energy per site:'), &
                       "The summary must label the per-site energy as such")
            call check(file_contains('test_t9_energy_summary.txt', '-2.000000000000'), &
                       "The per-site energy must be E_total / L = -8 / 4 = -2")
        end if

        call remove_file('test_t9_energy_summary.txt')
    end subroutine test_energy_per_site_uses_system_size


    !> Every output file must carry xc_smoothing_width and potential_tol
    !!
    !! A run with the V_xc discontinuity smoothed (w > 0) uses a MODIFIED XC
    !! functional: the energy per site of the reference case moves from -3.1938
    !! (w = 0.05) to -3.1719 (w = 0.2), i.e. 0.7%. Before this, a summary reading
    !! "SCF: CONVERGED / Final Energy: -3.193774" was indistinguishable between
    !! w = 0.05 and w = 0, so the archived numbers could not be told apart and the
    !! loss of parity with the C++ reference was invisible in the record.
    !!
    !! The test writes the same results twice, with w = 0 and w = 0.05, and
    !! requires that (a) both the tolerance and w appear in the summary and in the
    !! _density.dat / _convergence.dat headers, and (b) the w > 0 run carries an
    !! explicit statement that the functional was modified and has no C++ parity,
    !! which the w = 0 run must NOT carry.
    subroutine test_outputs_record_smoothing_and_tolerance()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use lsda_types, only: system_params_t
        use input_parser, only: input_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        type(input_params_t) :: inputs
        integer :: ierr

        sys_params%L = 3
        sys_params%Nup = 2
        sys_params%Ndown = 1
        sys_params%U = 4.0_dp
        sys_params%bc = 1

        results%converged = .true.
        results%n_iterations = 10
        results%final_density_error = 1.0e-10_dp
        results%final_potential_residual = 1.0e-9_dp
        results%final_energy = -5.0_dp

        allocate(results%density_up(3))
        allocate(results%density_down(3))
        allocate(results%eigvals(6))
        allocate(results%history%density_norms(2))
        allocate(results%history%energies(2))

        results%density_up = [0.7_dp, 0.7_dp, 0.6_dp]
        results%density_down = [0.3_dp, 0.3_dp, 0.4_dp]
        results%eigvals = [-1.0_dp, 0.0_dp, 1.0_dp, -0.5_dp, 0.5_dp, 1.5_dp]
        results%history%current_iter = 2
        results%history%density_norms = [1.0e-5_dp, 1.0e-10_dp]
        results%history%energies = [-4.8_dp, -5.0_dp]

        inputs%save_density = .true.
        inputs%save_eigenvalues = .false.
        inputs%store_history = .true.
        inputs%potential_tol = 1.0e-6_dp

        ! --- Unsmoothed run: exact C++ parity, no warning ---------------------
        inputs%output_prefix = 'test_prov_w0'
        inputs%xc_smoothing_width = 0.0_dp
        call write_results(results, sys_params, inputs, ierr)
        call check(ierr == ERROR_SUCCESS, "write_results (w = 0) should succeed")

        call check(file_contains('test_prov_w0_summary.txt', 'xc_smoothing_width'), &
                   "The summary must record xc_smoothing_width even when it is zero")
        call check(file_contains('test_prov_w0_summary.txt', 'potential_tol'), &
                   "The summary must record the potential tolerance actually used")
        call check(file_contains('test_prov_w0_density.dat', 'xc_smoothing_width'), &
                   "The density header must record xc_smoothing_width")
        call check(file_contains('test_prov_w0_convergence.dat', 'xc_smoothing_width'), &
                   "The convergence header must record xc_smoothing_width")
        call check(file_contains('test_prov_w0_density.dat', 'potential_tol'), &
                   "The density header must record the potential tolerance")
        call check(file_contains('test_prov_w0_convergence.dat', 'potential_tol'), &
                   "The convergence header must record the potential tolerance")
        call check(.not. file_contains('test_prov_w0_summary.txt', 'NO parity'), &
                   "An unsmoothed run must NOT be flagged as departing from the C++ reference")
        call check(.not. file_contains('test_prov_w0_summary.txt', 'not variational'), &
                   "At w = 0 V_xc IS the derivative of E_xc, so no variational caveat applies")

        ! --- Smoothed run: modified functional, must say so -------------------
        inputs%output_prefix = 'test_prov_w05'
        inputs%xc_smoothing_width = 0.05_dp
        call write_results(results, sys_params, inputs, ierr)
        call check(ierr == ERROR_SUCCESS, "write_results (w = 0.05) should succeed")

        call check(file_contains('test_prov_w05_summary.txt', '0.050000'), &
                   "The summary must record the actual value of w")
        call check(file_contains('test_prov_w05_summary.txt', 'MODIFIED'), &
                   "A smoothed run must state that the XC functional was modified")
        call check(file_contains('test_prov_w05_summary.txt', 'NO parity'), &
                   "A smoothed run must state that it has no parity with the C++ reference")
        call check(file_contains('test_prov_w05_density.dat', 'MODIFIED'), &
                   "The density header must carry the modified-functional warning")
        call check(file_contains('test_prov_w05_convergence.dat', 'MODIFIED'), &
                   "The convergence header must carry the modified-functional warning")

        ! The decisive line: E_xc is NOT smoothed. The smoothed V_xc reaches
        ! compute_total_energy only through the diagonalised V_eff, while E_xc
        ! is the unsmoothed get_exc, so with w > 0 the pair
        ! stops being (functional, derivative): the density the cycle converges
        ! to, n_w, is the stationary point of the SMOOTHED problem, not of the
        ! functional being evaluated. The reported energy is nevertheless a clean
        ! evaluation of the unsmoothed functional at n_w, and n_0 (the unsmoothed
        ! minimiser) IS a stationary point of it WHERE THAT FUNCTIONAL IS
        ! DIFFERENTIABLE, so there the error is SECOND order:
        ! E[n_w] - E[n_0] = O(|n_w - n_0|^2). Exactly in the regime that motivates
        ! w > 0 that caveat bites: the unsmoothed E_xc has a kink at n = 1 (the
        ! V_xc discontinuity is its derivative jump), so for densities pinned at
        ! n = 1 the error is FIRST order in |n_w - n_0|, with coefficient given by
        ! the V_xc jump. Either way it is an uncontrolled error -
        ! nothing bounds |n_w - n_0| - which is why it must be stated. Saying
        ! only "the functional was MODIFIED" invites the reader to assume E_xc
        ! followed V_xc.
        call check(file_contains('test_prov_w05_summary.txt', 'E_xc is NOT smoothed'), &
                   "A smoothed run must state that E_xc was left unsmoothed")
        call check(file_contains('test_prov_w05_summary.txt', 'not variational'), &
                   "A smoothed run must state that the reported energy is not variational")
        call check(file_contains('test_prov_w05_density.dat', 'E_xc is NOT smoothed'), &
                   "The density header must carry the unsmoothed-E_xc statement")
        call check(file_contains('test_prov_w05_convergence.dat', 'E_xc is NOT smoothed'), &
                   "The convergence header must carry the unsmoothed-E_xc statement")

        call remove_file('test_prov_w0_summary.txt')
        call remove_file('test_prov_w0_density.dat')
        call remove_file('test_prov_w0_convergence.dat')
        call remove_file('test_prov_w05_summary.txt')
        call remove_file('test_prov_w05_density.dat')
        call remove_file('test_prov_w05_convergence.dat')

        deallocate(results%density_up)
        deallocate(results%density_down)
        deallocate(results%eigvals)
        deallocate(results%history%density_norms)
        deallocate(results%history%energies)
    end subroutine test_outputs_record_smoothing_and_tolerance


    !> Every output file must identify the external potential that produced it
    !!
    !! Before this, write_summary recorded L, Nup, Ndown, U, potential_tol and w
    !! and nothing else: a _density.dat produced with 50% random impurities was
    !! literally irreproducible from its own contents, because neither the
    !! potential type, nor its strength, nor its concentration, nor the random
    !! seed, nor the boundary condition, nor the mixing appeared anywhere. Since
    !! the reference case of this phase IS a disordered run, that was the worst
    !! remaining traceability hole.
    !!
    !! Two runs are written: one with a pinned seed and one with pot_seed = -1.
    !! Both must record the seed AS GIVEN (sentinel included) and, next to it,
    !! the `effective_seed` that was actually handed to the generator. The second
    !! one used to be flagged "NOT reproducible", which was true while the drawn
    !! seed was thrown away; now that app/main.f90 resolves the seed explicitly
    !! (resolve_random_seed) the file must instead tell the reader which integer
    !! to put in pot_seed to replay the realisation.
    subroutine test_outputs_record_external_potential_provenance()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use lsda_types, only: system_params_t
        use input_parser, only: input_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        type(input_params_t) :: inputs
        integer :: ierr

        sys_params%L = 3
        sys_params%Nup = 2
        sys_params%Ndown = 1
        sys_params%U = -4.0_dp
        sys_params%bc = 0

        results%converged = .true.
        results%n_iterations = 10
        results%final_density_error = 1.0e-10_dp
        results%final_potential_residual = 1.0e-9_dp
        results%final_energy = -5.0_dp

        allocate(results%density_up(3))
        allocate(results%density_down(3))
        allocate(results%eigvals(6))
        allocate(results%history%density_norms(2))
        allocate(results%history%energies(2))

        results%density_up = [0.7_dp, 0.7_dp, 0.6_dp]
        results%density_down = [0.3_dp, 0.3_dp, 0.4_dp]
        results%eigvals = [-1.0_dp, 0.0_dp, 1.0_dp, -0.5_dp, 0.5_dp, 1.5_dp]
        results%history%current_iter = 2
        results%history%density_norms = [1.0e-5_dp, 1.0e-10_dp]
        results%history%energies = [-4.8_dp, -5.0_dp]

        inputs%save_density = .true.
        inputs%save_eigenvalues = .true.
        inputs%store_history = .true.
        inputs%potential_tol = 1.0e-6_dp
        inputs%energy_tol = 1.0e-9_dp
        inputs%max_iter = 4321
        inputs%xc_smoothing_width = 0.0_dp
        inputs%potential_type = 'impurity'
        inputs%V0 = -4.0_dp
        inputs%concentration = 50.0_dp
        inputs%bc_type = 'open'
        inputs%phase = 0.0_dp
        inputs%mixing_alpha = 0.05_dp
        inputs%use_adaptive_mixing = .true.

        ! --- Reproducible run: the seed pins the disorder realisation ---------
        inputs%output_prefix = 'test_potprov_seed'
        inputs%pot_seed = 12345
        inputs%effective_seed = 12345
        call write_results(results, sys_params, inputs, ierr)
        call check(ierr == ERROR_SUCCESS, "write_results (pinned seed) should succeed")

        call check(file_contains('test_potprov_seed_summary.txt', 'potential_type = impurity'), &
                   "The summary must record the external potential type")
        call check(file_contains('test_potprov_seed_summary.txt', 'V0 ='), &
                   "The summary must record the potential strength V0")
        call check(file_contains('test_potprov_seed_summary.txt', 'concentration ='), &
                   "The summary must record the impurity concentration")
        call check(file_contains('test_potprov_seed_summary.txt', 'pot_seed = 12345'), &
                   "The summary must record the random seed that fixes the disorder")
        call check(file_contains('test_potprov_seed_summary.txt', 'bc = open'), &
                   "The summary must record the boundary condition")
        call check(file_contains('test_potprov_seed_summary.txt', 'mixing_alpha ='), &
                   "The summary must record the mixing weight")
        call check(file_contains('test_potprov_seed_summary.txt', 'use_adaptive_mixing = T'), &
                   "The summary must record whether the adaptive controller was on")

        ! The same block must travel with the data files, not only with the summary.
        call check(file_contains('test_potprov_seed_density.dat', 'pot_seed = 12345'), &
                   "The density header must record the random seed")
        call check(file_contains('test_potprov_seed_density.dat', 'potential_type = impurity'), &
                   "The density header must record the external potential type")
        call check(file_contains('test_potprov_seed_convergence.dat', 'pot_seed = 12345'), &
                   "The convergence header must record the random seed")

        ! _eigenvalues.dat used not to receive the block at all (write_eigenvalues
        ! did not even take `inputs`), so a spectrum was archived with no record
        ! of the potential or the tolerances that produced it.
        call check(file_contains('test_potprov_seed_eigenvalues.dat', 'pot_seed = 12345'), &
                   "The eigenvalue header must record the random seed")
        call check(file_contains('test_potprov_seed_eigenvalues.dat', 'potential_type = impurity'), &
                   "The eigenvalue header must record the external potential type")
        call check(file_contains('test_potprov_seed_eigenvalues.dat', 'potential_tol'), &
                   "The eigenvalue header must record the SCF tolerance")

        ! The rest of the SCF contract: the second convergence criterion, the
        ! iteration budget and the twist phase. Without energy_tol and max_iter
        ! a "NOT CONVERGED" file cannot be told apart from a run that was simply
        ! never given enough iterations.
        call check(file_contains('test_potprov_seed_summary.txt', 'energy_tol ='), &
                   "The summary must record the energy tolerance, the second criterion")
        call check(file_contains('test_potprov_seed_summary.txt', 'max_iter = 4321'), &
                   "The summary must record the iteration budget actually granted")
        call check(file_contains('test_potprov_seed_summary.txt', 'phase ='), &
                   "The summary must record the twist phase of the boundary condition")
        call check(file_contains('test_potprov_seed_density.dat', 'energy_tol ='), &
                   "The density header must record the energy tolerance")
        call check(file_contains('test_potprov_seed_density.dat', 'max_iter = 4321'), &
                   "The density header must record the iteration budget")

        call check(.not. file_contains('test_potprov_seed_summary.txt', 'NOT reproducible'), &
                   "A pinned seed must NOT be flagged as irreproducible")
        call check(file_contains('test_potprov_seed_summary.txt', 'effective_seed = 12345'), &
                   "A pinned seed is also the effective seed, and must be recorded as such")

        ! --- Clock-seeded run: the drawn seed is what identifies it -----------
        ! Regression for the irreproducibility hole: with pot_seed = -1 the file
        ! must now carry the integer that was actually used (effective_seed) and
        ! the instruction to feed it back. Before this, the block recorded only
        ! the sentinel -1, which identifies no realisation at all.
        inputs%output_prefix = 'test_potprov_noseed'
        inputs%pot_seed = -1
        inputs%effective_seed = 987654
        call write_results(results, sys_params, inputs, ierr)
        call check(ierr == ERROR_SUCCESS, "write_results (clock seed) should succeed")

        call check(file_contains('test_potprov_noseed_summary.txt', 'pot_seed = -1'), &
                   "pot_seed must be recorded exactly as given, sentinel included")
        call check(file_contains('test_potprov_noseed_summary.txt', 'effective_seed = 987654'), &
                   "The seed actually drawn must be recorded next to the sentinel")
        call check(file_contains('test_potprov_noseed_summary.txt', 'set pot_seed = 987654'), &
                   "The summary must say how to replay the clock-seeded realisation")
        call check(file_contains('test_potprov_noseed_density.dat', 'effective_seed = 987654'), &
                   "The density header must carry the effective seed too")
        call check(.not. file_contains('test_potprov_noseed_summary.txt', 'NOT reproducible'), &
                   "A clock-seeded run IS reproducible once the drawn seed is recorded")

        call remove_file('test_potprov_seed_summary.txt')
        call remove_file('test_potprov_seed_density.dat')
        call remove_file('test_potprov_seed_eigenvalues.dat')
        call remove_file('test_potprov_seed_convergence.dat')
        call remove_file('test_potprov_noseed_summary.txt')
        call remove_file('test_potprov_noseed_density.dat')
        call remove_file('test_potprov_noseed_eigenvalues.dat')
        call remove_file('test_potprov_noseed_convergence.dat')

        deallocate(results%density_up)
        deallocate(results%density_down)
        deallocate(results%eigvals)
        deallocate(results%history%density_norms)
        deallocate(results%history%energies)
    end subroutine test_outputs_record_external_potential_provenance


    !> AAH provenance must contain the three parameters actually used to form V(i)
    subroutine test_quasiperiodic_provenance_records_generator_params()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use lsda_types, only: system_params_t
        use input_parser, only: input_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        type(input_params_t) :: inputs
        integer :: ierr

        sys_params%L = 4
        sys_params%Nup = 2
        sys_params%Ndown = 2
        sys_params%U = 4.0_dp
        results%converged = .true.
        results%n_iterations = 1
        results%final_energy = -4.0_dp

        inputs%output_prefix = 'test_prov_aah'
        inputs%potential_type = 'quasiperiodic'
        inputs%aah_lambda = 2.5_dp
        inputs%aah_beta = 0.25_dp
        inputs%aah_phi = 1.25_dp
        ! A deliberately unrelated legacy field: it must not be claimed as an
        ! AAH parameter in the output.
        inputs%pot_width = 99.0_dp
        inputs%save_density = .false.
        inputs%save_eigenvalues = .false.
        inputs%store_history = .false.

        call write_results(results, sys_params, inputs, ierr)

        call check(ierr == ERROR_SUCCESS, "write_results for AAH should succeed")
        call check(file_contains('test_prov_aah_summary.txt', 'aah_lambda =     2.500000'), &
                   "AAH provenance must record lambda")
        call check(file_contains('test_prov_aah_summary.txt', 'aah_beta =     0.250000'), &
                   "AAH provenance must record beta")
        call check(file_contains('test_prov_aah_summary.txt', 'aah_phi =     1.250000'), &
                   "AAH provenance must record phi")
        call check(.not. file_contains('test_prov_aah_summary.txt', 'pot_width ='), &
                   "AAH provenance must not record irrelevant pot_width")

        call remove_file('test_prov_aah_summary.txt')
    end subroutine test_quasiperiodic_provenance_records_generator_params


    !> The provenance must record the parameters of the potential ACTUALLY used
    !!
    !! Recording only V0 and the concentration made whole families of runs
    !! untraceable: a harmonic trap is defined by its spring constant, a double
    !! barrier by its geometry (barrier width, well depth, well width) and a
    !! disordered run by its strength - none of which appeared anywhere. The
    !! block is type-directed rather than exhaustive, so this test also checks
    !! that irrelevant parameters stay OUT: a harmonic run must not carry a well
    !! depth, nor a seed that never fed any random number generator.
    subroutine test_provenance_records_type_specific_params()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use lsda_types, only: system_params_t
        use input_parser, only: input_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        type(input_params_t) :: inputs
        integer :: ierr

        sys_params%L = 3
        sys_params%Nup = 2
        sys_params%Ndown = 1
        sys_params%U = 4.0_dp
        sys_params%bc = 1

        results%converged = .true.
        results%n_iterations = 10
        results%final_energy = -5.0_dp

        inputs%save_density = .false.
        inputs%save_eigenvalues = .false.
        inputs%store_history = .false.

        ! --- Harmonic trap: the spring constant IS the potential --------------
        inputs%output_prefix = 'test_provtype_harm'
        inputs%potential_type = 'harmonic'
        inputs%spring_constant = 0.00125_dp
        inputs%pot_seed = -1
        call write_results(results, sys_params, inputs, ierr)
        call check(ierr == ERROR_SUCCESS, "write_results (harmonic) should succeed")

        call check(file_contains('test_provtype_harm_summary.txt', 'spring_constant ='), &
                   "A harmonic run must record the spring constant that defines it")
        call check(.not. file_contains('test_provtype_harm_summary.txt', 'well_depth'), &
                   "A harmonic run must not record a barrier well depth it never used")
        call check(.not. file_contains('test_provtype_harm_summary.txt', 'pot_seed'), &
                   "A deterministic potential must not record a seed that fed nothing")
        call check(.not. file_contains('test_provtype_harm_summary.txt', 'NOT reproducible'), &
                   "A deterministic potential is reproducible whatever pot_seed says")

        ! --- Double barrier: three geometric parameters beyond V0 -------------
        inputs%output_prefix = 'test_provtype_bar'
        inputs%potential_type = 'barrier_double'
        inputs%barrier_width = 3.0_dp
        inputs%well_depth = -3.0_dp
        inputs%well_width = 20.0_dp
        call write_results(results, sys_params, inputs, ierr)
        call check(ierr == ERROR_SUCCESS, "write_results (barrier_double) should succeed")

        call check(file_contains('test_provtype_bar_summary.txt', 'barrier_width ='), &
                   "A double barrier must record the barrier width")
        call check(file_contains('test_provtype_bar_summary.txt', 'well_depth ='), &
                   "A double barrier must record the well depth")
        call check(file_contains('test_provtype_bar_summary.txt', 'well_width ='), &
                   "A double barrier must record the well width")

        ! --- Uniform disorder: strength plus the seed that fixes the draw -----
        inputs%output_prefix = 'test_provtype_dis'
        inputs%potential_type = 'random_uniform'
        inputs%disorder_strength = 1.75_dp
        inputs%pot_seed = 777
        call write_results(results, sys_params, inputs, ierr)
        call check(ierr == ERROR_SUCCESS, "write_results (random_uniform) should succeed")

        call check(file_contains('test_provtype_dis_summary.txt', 'disorder_strength ='), &
                   "A disordered run must record the disorder strength")
        call check(file_contains('test_provtype_dis_summary.txt', 'pot_seed = 777'), &
                   "A disordered run must record the seed that fixes the realisation")

        call remove_file('test_provtype_harm_summary.txt')
        call remove_file('test_provtype_bar_summary.txt')
        call remove_file('test_provtype_dis_summary.txt')
    end subroutine test_provenance_records_type_specific_params


    !> _convergence.dat must say which system produced the series
    !!
    !! The convergence history used to start straight at potential_tol: L, Nup,
    !! Ndown and U appeared in the other three files (they were written by each
    !! writer from sys_params) but never in the history, so an archived residual
    !! series could not be attributed to a system - not even to a filling or a
    !! sign of U. The four identification lines now come from the single shared
    !! provenance block, so this test also requires that all four files carry
    !! the SAME lines: any one output file must be self-sufficient.
    subroutine test_convergence_header_identifies_system()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use lsda_types, only: system_params_t
        use input_parser, only: input_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        type(input_params_t) :: inputs
        integer :: ierr

        ! Values chosen so that no two of the four can be confused with one
        ! another, and so that the sign of U is visible.
        sys_params%L = 7
        sys_params%Nup = 5
        sys_params%Ndown = 2
        sys_params%U = -3.25_dp
        sys_params%bc = 0

        results%converged = .true.
        results%n_iterations = 11
        results%final_density_error = 1.0e-10_dp
        results%final_potential_residual = 1.0e-9_dp
        results%final_energy = -5.0_dp

        allocate(results%density_up(7))
        allocate(results%density_down(7))
        allocate(results%eigvals(14))
        allocate(results%history%density_norms(2))
        allocate(results%history%energies(2))

        results%density_up = 0.5_dp
        results%density_down = 0.25_dp
        results%eigvals = 0.0_dp
        results%history%current_iter = 2
        results%history%density_norms = [1.0e-5_dp, 1.0e-10_dp]
        results%history%energies = [-4.8_dp, -5.0_dp]

        inputs%output_prefix = 'test_convid'
        inputs%save_density = .true.
        inputs%save_eigenvalues = .true.
        inputs%store_history = .true.
        inputs%potential_type = 'uniform'

        call write_results(results, sys_params, inputs, ierr)
        call check(ierr == ERROR_SUCCESS, "write_results should succeed")

        call check(file_contains('test_convid_convergence.dat', '# L = 7'), &
                   "The convergence header must record the lattice size L")
        call check(file_contains('test_convid_convergence.dat', '# Nup = 5'), &
                   "The convergence header must record Nup")
        call check(file_contains('test_convid_convergence.dat', '# Ndown = 2'), &
                   "The convergence header must record Ndown")
        call check(file_contains('test_convid_convergence.dat', '# U = -3.2500'), &
                   "The convergence header must record U, sign included")

        ! The same four lines must still reach the other three files, i.e. the
        ! move into the shared block must not have dropped them anywhere.
        call check(file_contains('test_convid_summary.txt', 'L = 7') .and. &
                   file_contains('test_convid_summary.txt', 'Nup = 5') .and. &
                   file_contains('test_convid_summary.txt', 'Ndown = 2') .and. &
                   file_contains('test_convid_summary.txt', 'U = -3.2500'), &
                   "The summary must still identify the system")
        call check(file_contains('test_convid_density.dat', '# L = 7') .and. &
                   file_contains('test_convid_density.dat', '# Nup = 5') .and. &
                   file_contains('test_convid_density.dat', '# Ndown = 2') .and. &
                   file_contains('test_convid_density.dat', '# U = -3.2500'), &
                   "The density header must still identify the system")
        call check(file_contains('test_convid_eigenvalues.dat', '# L = 7') .and. &
                   file_contains('test_convid_eigenvalues.dat', '# Nup = 5') .and. &
                   file_contains('test_convid_eigenvalues.dat', '# Ndown = 2') .and. &
                   file_contains('test_convid_eigenvalues.dat', '# U = -3.2500'), &
                   "The eigenvalue header must still identify the system")

        call remove_file('test_convid_summary.txt')
        call remove_file('test_convid_density.dat')
        call remove_file('test_convid_eigenvalues.dat')
        call remove_file('test_convid_convergence.dat')

        deallocate(results%density_up)
        deallocate(results%density_down)
        deallocate(results%eigvals)
        deallocate(results%history%density_norms)
        deallocate(results%history%energies)
    end subroutine test_convergence_header_identifies_system


    !> A random_uniform run must not claim a gaussian distribution
    !!
    !! The provenance used to write `distribution` for both random potentials,
    !! but the generator is chosen exclusively by potential_type
    !! (potential_factory never reads `distribution`). Since that field defaults
    !! to 'gaussian', a run declared as
    !!
    !!     potential_type = 'random_uniform'
    !!
    !! was correctly generated from the uniform generator while every one of its
    !! output files stated "distribution = gaussian": an output file asserting
    !! something objectively false about the calculation that produced it.
    !!
    !! The field is no longer written at all, because potential_type already
    !! identifies the generator completely; the `distribution` input key has
    !! since been removed from the namelist as well. The test requires that the
    !! word "gaussian" appear NOWHERE in any of the four files of a
    !! uniform-disorder run, while the generator itself remains identified.
    subroutine test_random_uniform_no_false_distribution()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use lsda_types, only: system_params_t
        use input_parser, only: input_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        type(input_params_t) :: inputs
        integer :: ierr

        sys_params%L = 3
        sys_params%Nup = 2
        sys_params%Ndown = 1
        sys_params%U = 4.0_dp
        sys_params%bc = 0

        results%converged = .true.
        results%n_iterations = 10
        results%final_density_error = 1.0e-10_dp
        results%final_potential_residual = 1.0e-9_dp
        results%final_energy = -5.0_dp

        allocate(results%density_up(3))
        allocate(results%density_down(3))
        allocate(results%eigvals(6))
        allocate(results%history%density_norms(2))
        allocate(results%history%energies(2))

        results%density_up = [0.7_dp, 0.7_dp, 0.6_dp]
        results%density_down = [0.3_dp, 0.3_dp, 0.4_dp]
        results%eigvals = [-1.0_dp, 0.0_dp, 1.0_dp, -0.5_dp, 0.5_dp, 1.5_dp]
        results%history%current_iter = 2
        results%history%density_norms = [1.0e-5_dp, 1.0e-10_dp]
        results%history%energies = [-4.8_dp, -5.0_dp]

        inputs%output_prefix = 'test_unifdist'
        inputs%save_density = .true.
        inputs%save_eigenvalues = .true.
        inputs%store_history = .true.
        inputs%potential_type = 'random_uniform'
        inputs%disorder_strength = 1.75_dp
        inputs%pot_seed = 777
        ! There is no `distribution` field to set any more; its old default,
        ! 'gaussian', is exactly what used to be written here.

        call write_results(results, sys_params, inputs, ierr)
        call check(ierr == ERROR_SUCCESS, "write_results (random_uniform) should succeed")

        call check(.not. file_contains('test_unifdist_summary.txt', 'gaussian'), &
                   "A uniform-disorder summary must not mention a gaussian distribution")
        call check(.not. file_contains('test_unifdist_density.dat', 'gaussian'), &
                   "A uniform-disorder density header must not mention gaussian")
        call check(.not. file_contains('test_unifdist_eigenvalues.dat', 'gaussian'), &
                   "A uniform-disorder eigenvalue header must not mention gaussian")
        call check(.not. file_contains('test_unifdist_convergence.dat', 'gaussian'), &
                   "A uniform-disorder convergence header must not mention gaussian")

        call check(.not. file_contains('test_unifdist_summary.txt', 'distribution'), &
                   "The unused `distribution` input must not be reported as provenance")

        ! What is written instead must still pin the generator and the draw.
        call check(file_contains('test_unifdist_summary.txt', &
                                 'potential_type = random_uniform'), &
                   "potential_type must identify the disorder generator")
        call check(file_contains('test_unifdist_summary.txt', 'disorder_strength ='), &
                   "The disorder strength must still be recorded")
        call check(file_contains('test_unifdist_summary.txt', 'pot_seed = 777'), &
                   "The seed fixing the realisation must still be recorded")

        call remove_file('test_unifdist_summary.txt')
        call remove_file('test_unifdist_density.dat')
        call remove_file('test_unifdist_eigenvalues.dat')
        call remove_file('test_unifdist_convergence.dat')

        deallocate(results%density_up)
        deallocate(results%density_down)
        deallocate(results%eigvals)
        deallocate(results%history%density_norms)
        deallocate(results%history%energies)
    end subroutine test_random_uniform_no_false_distribution


    !> Test density profile file format
    subroutine test_density_profile_format()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use lsda_types, only: system_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        integer :: ierr, io_unit, io_stat, site_read
        real(dp) :: nup_read, ndn_read, ntot_read
        character(len=256) :: line

        sys_params%L = 2
        sys_params%Nup = 1
        sys_params%Ndown = 1
        sys_params%U = 4.0_dp

        allocate(results%density_up(2))
        allocate(results%density_down(2))

        results%density_up = [0.5_dp, 0.5_dp]
        results%density_down = [0.5_dp, 0.5_dp]

        call write_density_profile(results, sys_params, 'test_fmt', ierr)

        ! Read file and check format
        open(newunit=io_unit, file='test_fmt_density.dat', status='old', iostat=io_stat)
        call check(io_stat == 0, "File should open successfully")

        ! Skip header lines (start with #)
        do
            read(io_unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) exit
            if (line(1:1) /= '#') then
                backspace(io_unit)
                exit
            end if
        end do

        ! Read first data line
        read(io_unit, *, iostat=io_stat) site_read, nup_read, ndn_read, ntot_read

        call check(io_stat == 0, "Data line should read successfully")
        call check(site_read == 1, "First site should be 1")
        call check(abs(nup_read - 0.5_dp) < 1.0e-6_dp, "n_up should match")
        call check(abs(ndn_read - 0.5_dp) < 1.0e-6_dp, "n_down should match")
        call check(abs(ntot_read - 1.0_dp) < 1.0e-6_dp, "n_total should be sum")

        close(io_unit, status='delete')

        deallocate(results%density_up)
        deallocate(results%density_down)
    end subroutine test_density_profile_format


    !> Test eigenvalues file format for occupied states
    subroutine test_eigenvalues_format_occupied()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use lsda_types, only: system_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        integer :: ierr, io_unit, io_stat, idx_read
        real(dp) :: eig_read
        character(len=256) :: line, spin_read, occ_read

        sys_params%L = 2
        sys_params%Nup = 1
        sys_params%Ndown = 1
        sys_params%U = 4.0_dp

        allocate(results%eigvals(4))
        results%eigvals = [-1.0_dp, 1.0_dp, -0.5_dp, 0.5_dp]

        call write_eigenvalues(results, sys_params, 'test_eig_occ', ierr)

        ! Read file and check format
        open(newunit=io_unit, file='test_eig_occ_eigenvalues.dat', status='old', iostat=io_stat)
        call check(io_stat == 0, "File should open successfully")

        ! Skip header
        do
            read(io_unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) exit
            if (line(1:1) /= '#') then
                backspace(io_unit)
                exit
            end if
        end do

        ! Read first eigenvalue (should be spin-up, occupied)
        read(io_unit, *, iostat=io_stat) idx_read, spin_read, eig_read, occ_read

        call check(io_stat == 0, "First eigenvalue should read")
        call check(idx_read == 1, "Index should be 1")
        call check(trim(adjustl(spin_read)) == 'up', "Spin should be 'up'")
        call check(abs(eig_read - (-1.0_dp)) < 1.0e-6_dp, "Eigenvalue should match")
        call check(trim(adjustl(occ_read)) == 'yes', "Should be occupied")

        close(io_unit, status='delete')

        deallocate(results%eigvals)
    end subroutine test_eigenvalues_format_occupied


    !> Test eigenvalues file format for unoccupied states
    subroutine test_eigenvalues_format_unoccupied()
        use fortuno_serial, only: check => serial_check
        use output_writer
        use lsda_types, only: system_params_t
        use kohn_sham_cycle, only: scf_results_t
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        type(scf_results_t) :: results
        type(system_params_t) :: sys_params
        integer :: ierr, io_unit, io_stat, idx_read
        real(dp) :: eig_read
        character(len=256) :: line, spin_read, occ_read

        sys_params%L = 3
        sys_params%Nup = 1  ! Only 1 occupied
        sys_params%Ndown = 1
        sys_params%U = 4.0_dp

        allocate(results%eigvals(6))
        results%eigvals = [-1.0_dp, 0.0_dp, 1.0_dp, -0.5_dp, 0.5_dp, 1.5_dp]

        call write_eigenvalues(results, sys_params, 'test_eig_unocc', ierr)

        ! Read file
        open(newunit=io_unit, file='test_eig_unocc_eigenvalues.dat', status='old', iostat=io_stat)
        call check(io_stat == 0, "File should open")

        ! Skip header
        do
            read(io_unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) exit
            if (line(1:1) /= '#') then
                backspace(io_unit)
                exit
            end if
        end do

        ! Read first eigenvalue (occupied)
        read(io_unit, *, iostat=io_stat) idx_read, spin_read, eig_read, occ_read
        call check(trim(adjustl(occ_read)) == 'yes', "First should be occupied")

        ! Read second eigenvalue (unoccupied)
        read(io_unit, *, iostat=io_stat) idx_read, spin_read, eig_read, occ_read
        call check(trim(adjustl(occ_read)) == 'no', "Second should be unoccupied")

        close(io_unit, status='delete')

        deallocate(results%eigvals)
    end subroutine test_eigenvalues_format_unoccupied

end program test_output_writer
