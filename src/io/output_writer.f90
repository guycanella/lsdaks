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
    use boundary_conditions, only: BC_TWISTED
    use lsda_errors, only: ERROR_SUCCESS, ERROR_FILE_WRITE
    implicit none
    private

    public :: write_results
    public :: write_density_profile
    public :: write_eigenvalues
    public :: write_convergence_history
    public :: write_summary

contains

    !> Write the provenance block shared by every output file
    !!
    !! The block records the inputs that define the calculation, so that a
    !! number found in an archived file can be traced back to the run that
    !! produced it. It records five groups:
    !!
    !!   * the system itself: L, Nup, Ndown and U. These are written here, from
    !!     `sys_params`, rather than by each caller, so that every output file
    !!     carries the SAME identification block and any one of them is
    !!     self-sufficient. _convergence.dat used to start straight at
    !!     potential_tol, i.e. an archived convergence series did not say which
    !!     system produced it;
    !!   * the SCF contract: potential_tol and energy_tol (the two convergence
    !!     criteria the run actually converged - or failed - to), max_iter,
    !!     mixing_alpha and use_adaptive_mixing;
    !!   * the XC functional: xc_smoothing_width (w). w > 0 is NOT a controlled
    !!     perturbation: it modifies the functional, and the energy per site
    !!     moves by ~0.7% between w = 0.05 and w = 0.2. Without this block a
    !!     summary reading "CONVERGED / Final Energy: -3.193774" would be
    !!     indistinguishable between a smoothed and an unsmoothed run, so an
    !!     explicit loss-of-parity warning is emitted whenever w > 0. The
    !!     warning also states that E_xc is NOT smoothed, i.e. that with w > 0
    !!     V_xc is no longer the functional derivative of the E_xc entering the
    !!     reported total energy (compute_total_energy calls the smoothed
    !!     `get_vxc` and the unsmoothed `get_exc`), so the reported energy is
    !!     not variational. The energy reported is a clean evaluation of the
    !!     UNSMOOTHED functional at the converged density n_w. WHERE THAT
    !!     FUNCTIONAL IS DIFFERENTIABLE, n_0 is a stationary point of it and
    !!     the error is SECOND order: E[n_w] - E[n_0] = O(|n_w - n_0|^2). That
    !!     is NOT the whole story in the regime that motivates w > 0: the
    !!     unsmoothed E_xc has a kink at n = 1 (the V_xc discontinuity IS its
    !!     derivative jump), so for sites pinned at n = 1 the error is FIRST
    !!     order in |n_w - n_0|, with coefficient given by the V_xc jump;
    !!   * the external potential: potential_type plus the parameters that
    !!     actually shape it for that type (V0, spring_constant, barrier and
    !!     well geometry, disorder strength, impurity positions,
    !!     concentration), and the seed for the stochastic types. A disordered
    !!     run whose seed is not recorded is irreproducible from its own file;
    !!   * the lattice boundary: bc and the twist phase.
    !!
    !! THE SEED. Two values are written for the stochastic potentials:
    !! `pot_seed` exactly as it was given (including the sentinel -1, "draw one
    !! from the system clock") and `effective_seed`, the non-negative integer
    !! that was actually handed to the generator. The second one is what
    !! identifies the disorder realisation: `app/main.f90` resolves the seed
    !! through `resolve_random_seed` before calling the generator, so setting
    !! `pot_seed = <effective_seed>` in the input replays the run bit for bit.
    !! Previously only `pot_seed` was recorded and a clock-seeded run was
    !! irreproducible even for the person who ran it.
    !!
    !! @param[in] io_unit    Open output unit
    !! @param[in] comment    Comment marker for the file format ("# " or "")
    !! @param[in] inputs     Full input record of the run (optional; when absent
    !!                       no provenance beyond the system identification is
    !!                       written)
    !! @param[in] sys_params System actually solved (optional; when absent the
    !!                       system identification lines are omitted)
    subroutine write_provenance(io_unit, comment, inputs, sys_params)
        integer, intent(in) :: io_unit
        character(len=*), intent(in) :: comment
        type(input_params_t), intent(in), optional :: inputs
        type(system_params_t), intent(in), optional :: sys_params

        logical :: seed_matters

        ! --- System identification ----------------------------------------
        ! Taken from sys_params (what was actually solved), not from the raw
        ! input record, so the block cannot drift from the calculation.
        if (present(sys_params)) then
            write(io_unit, '(A,I0)') comment // "L = ", sys_params%L
            write(io_unit, '(A,I0)') comment // "Nup = ", sys_params%Nup
            write(io_unit, '(A,I0)') comment // "Ndown = ", sys_params%Ndown
            write(io_unit, '(A,F0.4)') comment // "U = ", sys_params%U
        end if

        if (.not. present(inputs)) return

        ! --- SCF contract -------------------------------------------------
        write(io_unit, '(A,ES12.4)') comment // "potential_tol = ", inputs%potential_tol
        write(io_unit, '(A,ES12.4)') comment // "energy_tol = ", inputs%energy_tol
        write(io_unit, '(A,I0)') comment // "max_iter = ", inputs%max_iter
        write(io_unit, '(A,F10.6)') comment // "mixing_alpha = ", inputs%mixing_alpha
        write(io_unit, '(A,L1)') comment // "use_adaptive_mixing = ", inputs%use_adaptive_mixing

        ! --- External potential and boundary condition --------------------
        !
        ! Only the parameters that the requested potential_type actually reads
        ! are written: dumping every field indiscriminately would record, say, a
        ! well depth for a harmonic trap and invite the reader to believe it
        ! meant something.
        write(io_unit, '(A,A)') comment // "potential_type = ", trim(inputs%potential_type)
        write(io_unit, '(A,F12.6)') comment // "V0 = ", inputs%V0

        seed_matters = .false.

        select case (trim(inputs%potential_type))
        case ('harmonic')
            write(io_unit, '(A,F12.6)') comment // "spring_constant = ", inputs%spring_constant
        case ('random_uniform', 'random_gaussian')
            ! No `distribution` line. The generator is selected exclusively by
            ! potential_type (random_uniform vs random_gaussian), which
            ! identifies it completely. There used to be a `distribution` input
            ! key that nothing ever read; recording it stated a distribution the
            ! run did not use - its default was 'gaussian', so a random_uniform
            ! run was correctly generated as uniform while its own provenance
            ! claimed gaussian. The key has since been removed from the namelist.
            write(io_unit, '(A,F12.6)') comment // "disorder_strength = ", inputs%disorder_strength
            seed_matters = .true.
        case ('barrier_single')
            write(io_unit, '(A,I0)') comment // "position = ", inputs%position
            write(io_unit, '(A,I0)') comment // "width = ", inputs%width
        case ('barrier_double')
            write(io_unit, '(A,F12.6)') comment // "barrier_width = ", inputs%barrier_width
            write(io_unit, '(A,F12.6)') comment // "well_depth = ", inputs%well_depth
            write(io_unit, '(A,F12.6)') comment // "well_width = ", inputs%well_width
        case ('impurity')
            write(io_unit, '(A,F10.4)') comment // "concentration = ", inputs%concentration
            seed_matters = .true.
        case ('impurity_single')
            write(io_unit, '(A,F12.6)') comment // "pot_center = ", inputs%pot_center
        case ('impurity_multiple')
            write(io_unit, '(A,A)') comment // "imp_positions = ", trim(inputs%imp_positions_str)
        case ('quasiperiodic')
            ! These are the three values consumed by
            ! apply_potential_quasiperiodic.  `pot_width` belongs to no part of
            ! the AAH construction, so recording it here would make the output
            ! look reproducible while omitting the actual potential.
            write(io_unit, '(A,F12.6)') comment // "aah_lambda = ", inputs%aah_lambda
            write(io_unit, '(A,F12.6)') comment // "aah_beta = ", inputs%aah_beta
            write(io_unit, '(A,F12.6)') comment // "aah_phi = ", inputs%aah_phi
        end select

        ! The seed only identifies anything for the stochastic potentials; for a
        ! deterministic one it is noise in the record.
        if (seed_matters) then
            write(io_unit, '(A,I0)') comment // "pot_seed = ", inputs%pot_seed
            write(io_unit, '(A,I0)') comment // "effective_seed = ", inputs%effective_seed
            if (inputs%pot_seed < 0) then
                write(io_unit, '(A)') comment // &
                    "NOTE: pot_seed < 0 means the disorder realisation was drawn from the"
                write(io_unit, '(A)') comment // &
                    "      system clock. The seed actually used is effective_seed above:"
                write(io_unit, '(A,I0,A)') comment // &
                    "      set pot_seed = ", inputs%effective_seed, " to reproduce this run."
            end if
        end if

        write(io_unit, '(A,A)') comment // "bc = ", trim(inputs%bc_type)
        ! The unit is spelled out: `inputs%phase` is the value as typed by the
        ! user, in units of pi, while the twist angle that enters the
        ! Hamiltonian is phase*pi radians. Without the annotation a reader of
        ! this header cannot reproduce the run.
        write(io_unit, '(A,F12.6)') comment // "phase = ", inputs%phase
        write(io_unit, '(A)') comment // "  (phase is in units of pi: theta = phase*pi radians)"

        ! --- XC functional -------------------------------------------------
        ! F10.6 (not F0.6) so that the value keeps its leading zero: "0.050000"
        write(io_unit, '(A,F10.6)') comment // "xc_smoothing_width = ", inputs%xc_smoothing_width
        if (inputs%xc_smoothing_width > 0.0_dp) then
            write(io_unit, '(A)') comment // &
                "WARNING: the XC functional was MODIFIED (V_xc discontinuity at n = 1"
            write(io_unit, '(A)') comment // &
                "         linearly smoothed over the window [1-w, 1+w])."
            write(io_unit, '(A)') comment // &
                "         These results have NO parity with the C++ reference (w = 0)."
            write(io_unit, '(A)') comment // &
                "         E_xc is NOT smoothed: with w > 0, V_xc is not the functional"
            write(io_unit, '(A)') comment // &
                "         derivative of the E_xc used in the total energy; the reported"
            write(io_unit, '(A)') comment // &
                "         energy is not variational. Its error in the density deviation"
            write(io_unit, '(A)') comment // &
                "         from the unsmoothed minimiser is second order where E_xc is"
            write(io_unit, '(A)') comment // &
                "         differentiable, but first order in |dn| for densities pinned"
            write(io_unit, '(A)') comment // &
                "         at n = 1, with coefficient given by the V_xc jump."
        end if
    end subroutine write_provenance

    !> Write all simulation results
    !!
    !! Every requested file is attempted, even after one of them fails. The
    !! routine used to bail out at the first error, so a summary file that could
    !! not be opened (read-only directory, full disk, bad prefix) also threw away
    !! the density profile, the eigenvalues and the convergence history of a
    !! calculation that had already been performed. The first error encountered
    !! is the one reported: the caller only needs to know that the record on disk
    !! is incomplete, and it must know that whether the failure was the first or
    !! the last file.
    !!
    !! Nothing here depends on `results%converged`. A non-convergent run gets the
    !! same set of files as a convergent one - it is precisely the run whose
    !! density and spectrum need inspecting - and the files themselves state the
    !! convergence status (see write_provenance).
    !!
    !! @param[in] results      SCF results
    !! @param[in] sys_params   System parameters
    !! @param[in] inputs       Input parameters
    !! @param[out] ierr        Error code (first failure, if any)
    subroutine write_results(results, sys_params, inputs, ierr)
        type(scf_results_t), intent(in) :: results
        type(system_params_t), intent(in) :: sys_params
        type(input_params_t), intent(in) :: inputs
        integer, intent(out) :: ierr

        character(len=256) :: prefix
        integer :: ierr_step

        ierr = ERROR_SUCCESS
        prefix = trim(inputs%output_prefix)

        call write_summary(results, sys_params, inputs, ierr_step)
        if (ierr_step /= ERROR_SUCCESS .and. ierr == ERROR_SUCCESS) ierr = ierr_step

        if (inputs%save_density) then
            call write_density_profile(results, sys_params, prefix, ierr_step, inputs=inputs)
            if (ierr_step /= ERROR_SUCCESS .and. ierr == ERROR_SUCCESS) ierr = ierr_step
        end if

        if (inputs%save_eigenvalues) then
            call write_eigenvalues(results, sys_params, prefix, ierr_step, inputs=inputs)
            if (ierr_step /= ERROR_SUCCESS .and. ierr == ERROR_SUCCESS) ierr = ierr_step
        end if

        if (inputs%store_history .and. allocated(results%history%density_norms)) then
            call write_convergence_history(results, prefix, ierr_step, inputs=inputs, &
                                           sys_params=sys_params)
            if (ierr_step /= ERROR_SUCCESS .and. ierr == ERROR_SUCCESS) ierr = ierr_step
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
        ! The twist phase is only meaningful under twisted BC. The literal 2 used
        ! here was BC_PERIODIC, so the phase was printed for periodic runs (where
        ! it is always 0) and hidden for the twisted ones that actually carry it.
        if (sys_params%bc == BC_TWISTED) then
            ! sys_params%phase is in RADIANS (the user input is in units of pi
            ! and is converted once, in convert_to_system_params); say so, since the value
            ! printed here used to carry no unit at all.
            print '(A,F0.6,A)', "  Phase:            ", sys_params%phase, " rad"
        end if
        print '(A)', ""
        
        print '(A)', "SCF Convergence:"
        if (results%converged) then
            print '(A)', "  Status:           ✓ CONVERGED"
        else
            print '(A)', "  Status:           ✗ NOT CONVERGED"
        end if
        print '(A,I0)', "  Iterations:       ", results%n_iterations
        print '(A,ES12.4)', "  Final |ΔV|:       ", results%final_potential_residual
        print '(A,ES12.4)', "  Final |Δn|:       ", results%final_density_error
        ! `final_energy` is the TOTAL energy. Both numbers are printed, each under
        ! its own label: the per-site value used to be obtained by dividing by
        ! size(results%density_up) and, when that array was not allocated, the
        ! TOTAL energy was printed verbatim under the label "per site" - off by a
        ! factor of L with no indication. The divisor is now sys_params%L, the
        ! system actually solved, which does not depend on which optional
        ! components of `results` the SCF cycle happened to fill in.
        print '(A,F20.12)', "  Final Total Energy:    ", results%final_energy
        if (sys_params%L > 0) then
            print '(A,F20.12)', "  Final Energy per site: ", &
                results%final_energy / real(sys_params%L, dp)
        end if
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
        call write_provenance(io_unit, "", inputs, sys_params)
        write(io_unit, '(A)') ""

        ! The human-readable line and a machine-readable flag. "SCF: CONVERGED"
        ! and "SCF: NOT CONVERGED" differ by a word in the middle, so any script
        ! grepping for the former also matches the latter unless it is careful;
        ! `converged = T|F` is unambiguous and is the field to test.
        if (results%converged) then
            write(io_unit, '(A)') "SCF: CONVERGED"
        else
            write(io_unit, '(A)') "SCF: NOT CONVERGED"
        end if
        write(io_unit, '(A,L1)') "converged = ", results%converged
        write(io_unit, '(A,I0)') "Iterations: ", results%n_iterations
        write(io_unit, '(A,ES12.4)') "Final |ΔV|: ", results%final_potential_residual
        write(io_unit, '(A,ES12.4)') "Final |Δn|: ", results%final_density_error
        ! Both energies, each labelled. "Final Energy:" was the energy per site
        ! without saying so, and fell back to the total energy - under the same
        ! label - whenever results%density_up was not allocated.
        write(io_unit, '(A,F20.12)') "Final Total Energy: ", results%final_energy
        if (sys_params%L > 0) then
            write(io_unit, '(A,F20.12)') "Final Energy per site: ", &
                results%final_energy / real(sys_params%L, dp)
        end if


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
    !! @param[in] inputs       Full input record of the run (optional). When
    !!                         given, the provenance block (SCF tolerances and
    !!                         mixing, external potential and seed, XC smoothing
    !!                         width) is written into the header, which is what
    !!                         makes a disordered run reproducible from the file.
    subroutine write_density_profile(results, sys_params, prefix, ierr, inputs)
        type(scf_results_t), intent(in) :: results
        type(system_params_t), intent(in) :: sys_params
        character(len=*), intent(in) :: prefix
        integer, intent(out) :: ierr
        type(input_params_t), intent(in), optional :: inputs

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
        ! `inputs` is optional here and optional there: passing an absent
        ! optional through is legal and simply means "no provenance block"
        ! beyond the system identification, which comes from sys_params.
        call write_provenance(io_unit, "# ", inputs, sys_params)
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
    !! @param[in] inputs       Full input record of the run (optional). When
    !!                         given, the provenance block is written into the
    !!                         header, exactly as in the density and convergence
    !!                         files: an eigenvalue spectrum is as meaningless
    !!                         without the potential and the tolerances that
    !!                         produced it as a density profile is.
    subroutine write_eigenvalues(results, sys_params, prefix, ierr, inputs)
        type(scf_results_t), intent(in) :: results
        type(system_params_t), intent(in) :: sys_params
        character(len=*), intent(in) :: prefix
        integer, intent(out) :: ierr
        type(input_params_t), intent(in), optional :: inputs

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
        call write_provenance(io_unit, "# ", inputs, sys_params)
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
    !! @param[in] results    SCF results
    !! @param[in] prefix     Output file prefix
    !! @param[out] ierr      Error code
    !! @param[in] inputs     Full input record of the run (optional). When given,
    !!                       the provenance block is written into the header.
    !! @param[in] sys_params System actually solved (optional). When given, the
    !!                       header identifies the system (L, Nup, Ndown, U)
    !!                       exactly as the other three output files do; without
    !!                       it the series cannot be attributed to a system.
    subroutine write_convergence_history(results, prefix, ierr, inputs, sys_params)
        type(scf_results_t), intent(in) :: results
        character(len=*), intent(in) :: prefix
        integer, intent(out) :: ierr
        type(input_params_t), intent(in), optional :: inputs
        type(system_params_t), intent(in), optional :: sys_params

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
        call write_provenance(io_unit, "# ", inputs, sys_params)
        write(io_unit, '(A)') "#"

        if (allocated(results%history%potential_residuals)) then
            write(io_unit, '(A)') "# |dV| is the potential self-consistency residual"
            write(io_unit, '(A)') "# (the actual convergence criterion; |dn| is diagnostic only)"
            write(io_unit, '(A)') "# Columns: iteration  |Δn|  energy  |ΔV|"

            do i = 1, results%history%current_iter
                write(io_unit, '(I6,3ES20.10)') i, &
                                                results%history%density_norms(i), &
                                                results%history%energies(i), &
                                                results%history%potential_residuals(i)
            end do
        else
            write(io_unit, '(A)') "# Columns: iteration  |Δn|  energy"

            do i = 1, results%history%current_iter
                write(io_unit, '(I6,2ES20.10)') i, &
                                                results%history%density_norms(i), &
                                                results%history%energies(i)
            end do
        end if
        
        close(io_unit)
        
        print '(A,A)', "  Convergence history written to: ", trim(filename)
        
        ierr = ERROR_SUCCESS
        
    end subroutine write_convergence_history

end module output_writer
