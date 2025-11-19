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
            test("write_eigenvalues_not_allocated", test_write_eigenvalues_not_allocated), &
            test("write_convergence_history", test_write_convergence_history), &
            test("write_convergence_not_allocated", test_write_convergence_not_allocated), &
            test("write_results_all_enabled", test_write_results_all_enabled), &
            test("write_results_density_disabled", test_write_results_density_disabled), &
            test("density_profile_format", test_density_profile_format), &
            test("eigenvalues_format_occupied", test_eigenvalues_format_occupied), &
            test("eigenvalues_format_unoccupied", test_eigenvalues_format_unoccupied) &
        ])
    end function get_output_writer_tests


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
