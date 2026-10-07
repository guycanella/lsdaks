!> Module for parsing input parameters from command line and namelist files
!!
!! Supports two input modes:
!! 1. Namelist file: fpm run lsdaks -- --input input.txt
!! 2. Command line:  fpm run lsdaks -- --L 10 --Nup 5 --Ndown 5 --U 4.0
!!
!! Priority: Command line arguments override namelist values
module input_parser
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
    use, intrinsic :: iso_fortran_env, only: IOSTAT_END, IOSTAT_EOR
    use lsda_constants, only: dp, PI, ITER_MAX, SCF_DENSITY_TOL, SCF_ENERGY_TOL, &
                              SCF_POTENTIAL_TOL, MIX_ALPHA
    use lsda_types, only: system_params_t
    use kohn_sham_cycle, only: scf_params_t
    use xc_lsda, only: XC_SMOOTHING_WIDTH_MAX
    use boundary_conditions, only: BC_OPEN, BC_PERIODIC, BC_TWISTED
    use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT, ERROR_FILE_NOT_FOUND
    implicit none
    private

    !> Default AAH modulation wavenumber: the inverse golden mean (√5-1)/2.
    !!
    !! Written as a literal instead of `(sqrt(5.0_dp) - 1.0_dp)/2.0_dp` because
    !! a default component initializer must be a constant expression.
    real(dp), parameter :: AAH_BETA_DEFAULT = 0.6180339887498948_dp

    !> Length of the read buffer used for `imp_positions_str`
    !!
    !! Several times the length of the field it is copied into, on purpose: it
    !! is what lets `read_namelist_file` SEE a value that does not fit instead
    !! of receiving it already truncated (see the declaration in that routine).
    !! A value longer than this buffer is truncated too, and is caught by the
    !! second guard there.
    integer, parameter :: IMP_STR_BUFFER_LEN = 4000

    !> Input parameters type
    type, public :: input_params_t
        ! System parameters
        integer :: L = 10
        integer :: Nup = 5
        integer :: Ndown = 5
        real(dp) :: U = 4.0_dp
        character(len=20) :: bc_type = 'periodic'
        !> Twist angle of the twisted boundary condition, **in units of π**.
        !!
        !! This is the user-facing unit and it matches the C++ reference, whose
        !! prompt reads "Give the phase of twisted boundary condition (in units
        !! of Pi)" and which multiplies the value by π right after reading it
        !! (`original/lsda_interface.cc:62-67,148-153`). The conversion to
        !! radians happens exactly once, in `convert_to_system_params`, when
        !! filling `system_params_t%phase`; everything downstream (the
        !! Hamiltonian, the BC validation) works in radians.
        real(dp) :: phase = 0.0_dp
        !> Directory holding the XC tables (`xc_table_u<|U|>.dat`).
        !!
        !! Relative paths are resolved against the working directory. If the
        !! requested table is not found here, `app/main.f90` also tries the
        !! directory named by the environment variable `LSDAKS_TABLE_DIR`, so
        !! the executable can be run from outside the repository root.
        character(len=256) :: table_dir = 'tables'

        ! Potential parameters
        character(len=20) :: potential_type = 'uniform'
        real(dp) :: V0 = 0.0_dp
        real(dp) :: pot_center = 0.0_dp
        real(dp) :: pot_width = 1.0_dp
        real(dp) :: spring_constant = 0.001_dp        ! For harmonic potential (k)
        real(dp) :: concentration = 50.0_dp           ! For random impurities (%)
        integer :: pot_seed = -1                      ! Random seed (-1 = system time)
        character(len=500) :: imp_positions_str = ''  ! For impurity_multiple: '10, 25, 40, 55'
        real(dp) :: disorder_strength = 2.0_dp        ! For random disorder (W or sigma)
        !> Seed actually used to draw the disorder realisation.
        !!
        !! NOT a namelist key: it is filled by `app/main.f90` after resolving
        !! `pot_seed` (see `resolve_random_seed`), so that a run started with
        !! `pot_seed = -1` still records the integer that identifies its
        !! realisation and can be replayed by feeding it back as `pot_seed`.
        integer :: effective_seed = -1
        ! Aubry-André-Harper quasiperiodic potential V(i) = lambda*cos(2*pi*beta*i + phi).
        ! There is no C++ counterpart for this potential, so the names and the
        ! defaults are ours: beta = (sqrt(5)-1)/2 is the inverse golden mean,
        ! the standard incommensurate modulation of the AAH model.
        real(dp) :: aah_lambda = 1.0_dp               ! Modulation amplitude (lambda)
        !> Modulation wavenumber β of the AAH potential, in units of 2π.
        !!
        !! @note β must be IRRATIONAL for the potential to be quasiperiodic.
        !! A simple rational value is accepted in silence and produces a
        !! genuinely PERIODIC potential of period 1/β sites - a different
        !! physical problem, with Bloch (extended) eigenstates instead of the
        !! localization transition the AAH model is used for: β = 0.5 gives a
        !! two-site alternation, β = 1/3 a three-site one, and no value of λ
        !! localizes them. The code cannot check this (no finite floating-point
        !! number is irrational, and every β is a rational with denominator
        !! 2^52), so it is the user's responsibility; the default is the inverse
        !! golden mean, the best-conditioned irrational for a finite lattice.
        !! A near-rational β (0.4999999) is just as suspect: the period is then
        !! longer than L and the modulation looks periodic over the lattice.
        real(dp) :: aah_beta = AAH_BETA_DEFAULT
        real(dp) :: aah_phi = 0.0_dp                  ! Modulation phase (phi), in radians
        integer :: position = 50                      ! For barrier_single: center position
        integer :: width = 5                          ! For barrier_single: width
        real(dp) :: barrier_width = 3.0_dp            ! For barrier_double: width of each barrier (Lb)
        real(dp) :: well_depth = -3.0_dp              ! For barrier_double: potential in well (Vwell, typically negative)
        real(dp) :: well_width = 20.0_dp              ! For barrier_double: width of well between barriers (Lwell)
        integer :: position1 = 35                     ! DEPRECATED for barrier_double
        integer :: width1 = 3                         ! DEPRECATED for barrier_double
        integer :: position2 = 65                     ! DEPRECATED for barrier_double
        integer :: width2 = 3                         ! DEPRECATED for barrier_double
        
        ! SCF parameters
        integer :: max_iter = ITER_MAX
        !> DIAGNOSTIC ONLY: ||Δn|| is reported but is no longer a convergence
        !! criterion (it scales with the mixing weight). Changing it does not
        !! tighten or loosen convergence; use potential_tol / energy_tol for that.
        real(dp) :: density_tol = SCF_DENSITY_TOL
        real(dp) :: energy_tol = SCF_ENERGY_TOL
        real(dp) :: potential_tol = SCF_POTENTIAL_TOL  ! Convergence tol for the potential residual
        real(dp) :: mixing_alpha = MIX_ALPHA
        logical :: verbose = .true.
        logical :: store_history = .true.
        logical :: use_adaptive_mixing = .true.   ! Use adaptive mixing (C++ behavior)
        !> Half-width w of the linear smoothing of V_xc around n = 1 (0 = off;
        !! must be < 1).
        !!
        !! POLICY (T20): the default is and stays 0, i.e. the exact BALDA
        !! functional of the C++ reference, including its discontinuity at
        !! n = 1. Smoothing is an explicit opt-in and is reported in the
        !! output header. Reason: w > 0 CHANGES THE FUNCTIONAL, not only the
        !! numerics - E/L moves by ~0.7% between w = 0.05 and w = 0.2 - so
        !! any self-consistent reference (the phase-2 harmonic-trap target in
        !! particular) must be produced with the original functional. Cases
        !! that only converge with w > 0 (a Mott plateau sitting on n = 1 in
        !! a trap or a double-barrier well) are documented per case, never
        !! made the default.
        real(dp) :: xc_smoothing_width = 0.0_dp

        ! Output parameters
        character(len=100) :: output_prefix = 'lsda_output'
        logical :: save_density = .true.
        logical :: save_eigenvalues = .true.
        logical :: save_wavefunction = .false.
    end type input_params_t

    public :: parse_inputs
    !> Public so that the namelist groups (and therefore every accepted key) can
    !! be exercised by the unit tests without spawning the executable.
    public :: read_namelist_file
    !> Public so that the gate that decides whether a namelist group is in the
    !! file can be compared, case by case, against what `read(nml=)` actually
    !! does with the same records. A divergence between the two is invisible in
    !! a run (a group the gate misses simply keeps its defaults), so it can only
    !! be caught by testing the two side by side.
    public :: namelist_group_present
    public :: validate_inputs
    public :: convert_to_system_params
    public :: convert_to_scf_params
    public :: parse_int_list
    public :: barrier_single_bounds

contains

    !> Parse inputs from command line or namelist file
    !!
    !! @param[out] inputs Parsed input parameters
    !! @param[out] ierr   Error code (0 = success)
    subroutine parse_inputs(inputs, ierr)
        type(input_params_t), intent(out) :: inputs
        integer, intent(out) :: ierr
        
        character(len=256) :: arg, input_file
        integer :: nargs, i
        logical :: use_input_file
        
        ierr = ERROR_SUCCESS
        use_input_file = .false.
        
        ! Initialize with defaults
        inputs = input_params_t()
        
        ! Check command line arguments
        nargs = command_argument_count()
        
        if (nargs == 0) then
            ! No arguments: try default input file
            input_file = 'input.txt'
            inquire(file=input_file, exist=use_input_file)
            if (.not. use_input_file) then
                print *, "WARNING: No input file or arguments provided."
                print *, "Using default parameters."
                return
            end if
        else
            ! Parse command line
            i = 1
            do while (i <= nargs)
                call get_command_argument(i, arg)
                
                select case (trim(arg))
                case ('--input', '-i')
                    if (i + 1 > nargs) then
                        print *, "ERROR: --input requires a filename"
                        ierr = ERROR_INVALID_INPUT
                        return
                    end if
                    call get_command_argument(i + 1, input_file)
                    use_input_file = .true.
                    i = i + 2
                    
                case ('--L')
                    call parse_integer_arg(i, nargs, inputs%L, ierr)
                    if (ierr /= ERROR_SUCCESS) return
                    i = i + 2
                    
                case ('--Nup')
                    call parse_integer_arg(i, nargs, inputs%Nup, ierr)
                    if (ierr /= ERROR_SUCCESS) return
                    i = i + 2
                    
                case ('--Ndown')
                    call parse_integer_arg(i, nargs, inputs%Ndown, ierr)
                    if (ierr /= ERROR_SUCCESS) return
                    i = i + 2
                    
                case ('--U')
                    call parse_real_arg(i, nargs, inputs%U, ierr)
                    if (ierr /= ERROR_SUCCESS) return
                    i = i + 2
                    
                case ('--bc')
                    if (i + 1 > nargs) then
                        print *, "ERROR: --bc requires a value (open/periodic/twisted)"
                        ierr = ERROR_INVALID_INPUT
                        return
                    end if
                    call get_command_argument(i + 1, inputs%bc_type)
                    i = i + 2
                    
                case ('--phase')
                    call parse_real_arg(i, nargs, inputs%phase, ierr)
                    if (ierr /= ERROR_SUCCESS) return
                    i = i + 2
                    
                case ('--potential')
                    if (i + 1 > nargs) then
                        print *, "ERROR: --potential requires a type"
                        ierr = ERROR_INVALID_INPUT
                        return
                    end if
                    call get_command_argument(i + 1, inputs%potential_type)
                    i = i + 2
                    
                case ('--V0')
                    call parse_real_arg(i, nargs, inputs%V0, ierr)
                    if (ierr /= ERROR_SUCCESS) return
                    i = i + 2
                    
                case ('--concentration')
                    call parse_real_arg(i, nargs, inputs%concentration, ierr)
                    if (ierr /= ERROR_SUCCESS) return
                    i = i + 2
                    
                case ('--verbose')
                    inputs%verbose = .true.
                    i = i + 1
                    
                case ('--quiet')
                    inputs%verbose = .false.
                    i = i + 1
                    
                case ('--help', '-h')
                    call print_help()
                    ierr = -1  ! Signal to exit gracefully
                    return
                    
                case default
                    print *, "WARNING: Unknown argument: ", trim(arg)
                    i = i + 1
                end select
            end do
        end if
        
        ! Read namelist file if specified
        if (use_input_file) then
            call read_namelist_file(input_file, inputs, ierr)
            if (ierr /= ERROR_SUCCESS) return
        end if
        
    end subroutine parse_inputs

    !> Read parameters from namelist file
    !!
    !! The file is first slurped into a character array (`read_input_records`)
    !! and the four groups are then read from that array as an INTERNAL file.
    !! This is deliberate: reading a namelist straight from the external unit
    !! makes the result depend on whether the file ends with a newline, because
    !! gfortran returns iostat = -1 ("End of file") for the group that closes an
    !! unterminated last line even though every value in it was assigned
    !! correctly. Since that status is indistinguishable from a genuinely broken
    !! group, a perfectly valid input file written by an editor or generator that
    !! omits the final newline was rejected as malformed. In an internal file
    !! every element is a complete record by construction, so the last group
    !! reads normally while a group with no closing '/' still hits end-of-file
    !! and is still rejected.
    !!
    !! @param[in]    filename Input file path
    !! @param[inout] inputs   Input parameters (updated from file)
    !! @param[out]   ierr     Error code
    subroutine read_namelist_file(filename, inputs, ierr)
        character(len=*), intent(in) :: filename
        type(input_params_t), intent(inout) :: inputs
        integer, intent(out) :: ierr

        integer :: io_stat
        character(len=512) :: io_msg
        character(len=:), allocatable :: records(:)
        logical :: file_exists

        ! System namelist variables
        integer :: L, Nup, Ndown
        real(dp) :: U, phase
        character(len=20) :: bc
        character(len=256) :: table_dir

        ! Potential namelist variables
        character(len=20) :: potential_type
        real(dp) :: V0, pot_center, pot_width, concentration, disorder_strength
        real(dp) :: barrier_width, well_depth, well_width
        real(dp) :: aah_lambda, aah_beta, aah_phi
        integer :: pot_seed, position, width, position1, width1, position2, width2
        !> Read buffer for the impurity list, DELIBERATELY much longer than the
        !! field it is copied into (`input_params_t%imp_positions_str`).
        !!
        !! gfortran truncates a namelist value that does not fit its variable
        !! SILENTLY (measured: iostat = 0, the tail simply dropped), and for a
        !! list of impurity sites that is the most expensive kind of input
        !! error: the run continues with fewer impurities than asked for and
        !! produces a perfectly plausible result for a different system.
        !!
        !! Reading into a buffer several times the size of the destination is
        !! what makes the overflow VISIBLE: the comparison below is
        !! `len_trim(buffer) > len(inputs%imp_positions_str)`, which needs the
        !! full value to be present. Testing "the destination field is full"
        !! instead does not work, because the cut can land on a blank (a list
        !! truncated after "..., 199, " leaves the 500th character blank).
        character(len=IMP_STR_BUFFER_LEN) :: imp_positions_str
        real(dp) :: spring_constant

        ! SCF namelist variables
        ! NOTE on /scf/: potential_tol and energy_tol are the convergence
        ! criteria; density_tol is accepted for backwards compatibility and
        ! reported, but it is DIAGNOSTIC ONLY and does not affect convergence.
        ! mixing_alpha is the weight of the new potential, and also the starting
        ! value of the adaptive controller when use_adaptive_mixing = .true.
        integer :: max_iter
        real(dp) :: density_tol, energy_tol, potential_tol, mixing_alpha
        real(dp) :: xc_smoothing_width
        logical :: verbose, store_history, use_adaptive_mixing

        ! Output namelist variables
        character(len=100) :: output_prefix
        logical :: save_density, save_eigenvalues, save_wavefunction

        namelist /system/ L, Nup, Ndown, U, bc, phase, table_dir
        ! NOTE on /potential/:
        !
        ! * `spring_constant` (the harmonic k) MUST be listed here. It was
        !   missing while the field existed in input_params_t and was consumed by
        !   app/main.f90, so `&potential ... spring_constant = 0.02 /` ran
        !   silently with the default k = 0.001 and made every comparison of the
        !   harmonic trap against the C++ reference meaningless. This is the
        !   namelist half of Bug #5 ("Harmonic Parameter Not Passed"), whose
        !   main.f90 half was fixed earlier.
        ! * `distribution` was REMOVED. Nothing ever read it: the generator is
        !   chosen by potential_type alone ('random_uniform' / 'random_gaussian').
        !   An input file that still sets it is now rejected by
        !   check_namelist_status, which prints a dedicated migration hint.
        namelist /potential/ potential_type, V0, pot_center, pot_width, &
                             concentration, pot_seed, imp_positions_str, &
                             disorder_strength, spring_constant, position, width, &
                             barrier_width, well_depth, well_width, &
                             aah_lambda, aah_beta, aah_phi, &
                             position1, width1, position2, width2

        namelist /scf/ max_iter, density_tol, energy_tol, potential_tol, mixing_alpha, &
                       verbose, store_history, use_adaptive_mixing, xc_smoothing_width

        namelist /output/ output_prefix, save_density, save_eigenvalues, &
                                                            save_wavefunction
        
        ierr = ERROR_SUCCESS
        
        ! Check if file exists
        inquire(file=filename, exist=file_exists)
        if (.not. file_exists) then
            print *, "ERROR: Input file not found: ", trim(filename)
            ierr = ERROR_FILE_NOT_FOUND
            return
        end if
        
        ! Initialize namelist variables with current values
        L = inputs%L
        Nup = inputs%Nup
        Ndown = inputs%Ndown
        U = inputs%U
        bc = inputs%bc_type
        phase = inputs%phase
        table_dir = inputs%table_dir

        potential_type = inputs%potential_type
        V0 = inputs%V0
        pot_center = inputs%pot_center
        pot_width = inputs%pot_width
        concentration = inputs%concentration
        pot_seed = inputs%pot_seed
        imp_positions_str = inputs%imp_positions_str
        disorder_strength = inputs%disorder_strength
        spring_constant = inputs%spring_constant
        aah_lambda = inputs%aah_lambda
        aah_beta = inputs%aah_beta
        aah_phi = inputs%aah_phi
        position = inputs%position
        width = inputs%width
        barrier_width = inputs%barrier_width
        well_depth = inputs%well_depth
        well_width = inputs%well_width
        position1 = inputs%position1
        width1 = inputs%width1
        position2 = inputs%position2
        width2 = inputs%width2
        
        max_iter = inputs%max_iter
        density_tol = inputs%density_tol
        energy_tol = inputs%energy_tol
        potential_tol = inputs%potential_tol
        mixing_alpha = inputs%mixing_alpha
        verbose = inputs%verbose
        store_history = inputs%store_history
        use_adaptive_mixing = inputs%use_adaptive_mixing
        xc_smoothing_width = inputs%xc_smoothing_width
        
        output_prefix = inputs%output_prefix
        save_density = inputs%save_density
        save_eigenvalues = inputs%save_eigenvalues
        save_wavefunction = inputs%save_wavefunction
        
        ! Slurp the file into complete records (see the note on this routine).
        call read_input_records(filename, records, ierr)
        if (ierr /= ERROR_SUCCESS) return

        ! Read the four groups. A MISSING group is legitimate (each group is
        ! optional and its keys keep their defaults) and only warns; anything
        ! else aborts. Ignoring iostat, as this routine used to do, turned every
        ! misspelled key and every malformed value into a silent fallback to the
        ! default - the most expensive kind of input error, because the run
        ! completes and produces plausible numbers for a system nobody asked for.
        !
        ! Presence of the group is decided BEFORE the read, by scanning the
        ! records for the header, and only a group that is really there is read.
        ! That keeps the "absent group" note independent of how the compiler
        ! reports an absent group (gfortran returns iostat = 0 on an internal
        ! file and iostat = -1 on an external one), and it leaves a single
        ! meaning for a non-zero iostat: the group is in the file and is broken.
        !
        ! No `rewind` between the reads: every read of an internal file starts at
        ! its first record.
        if (namelist_group_present(records, 'system')) then
            io_msg = ''
            read(records, nml=system, iostat=io_stat, iomsg=io_msg)
            call check_namelist_status('system', filename, io_stat, io_msg, ierr)
            if (ierr /= ERROR_SUCCESS) return
        else
            call print_missing_group_note('system', filename)
        end if

        if (namelist_group_present(records, 'potential')) then
            io_msg = ''
            read(records, nml=potential, iostat=io_stat, iomsg=io_msg)
            call check_namelist_status('potential', filename, io_stat, io_msg, ierr)
            if (ierr /= ERROR_SUCCESS) return
            ! An impurity list that does not fit its destination field must be
            ! REJECTED, not quietly shortened. gfortran truncates an oversized
            ! namelist value in silence (iostat = 0, tail dropped), and a
            ! truncated list of sites still yields a perfectly plausible - and
            ! wrong - calculation. parse_int_list rejects every other malformed
            ! list (empty field, non-integer token, out of range, duplicate), so
            ! tolerating the loss of a tail by field size was incoherent.
            !
            ! Two guards, because the value was read into a buffer several times
            ! longer than the field (see its declaration):
            !   1. the buffer itself overflowed - only for a truly absurd list;
            !   2. the value fits the buffer but not the field.
            if (len_trim(imp_positions_str) >= len(imp_positions_str)) then
                print *, "ERROR reading &potential in " // trim(filename) // ":"
                print '(A,I0,A)', "  imp_positions_str is longer than the ", &
                                  len(imp_positions_str), "-character read buffer and was truncated."
                ierr = ERROR_INVALID_INPUT
                return
            end if

            if (len_trim(imp_positions_str) > len(inputs%imp_positions_str)) then
                print *, "ERROR reading &potential in " // trim(filename) // ":"
                print '(A,I0,A,I0,A)', "  imp_positions_str has ", len_trim(imp_positions_str), &
                                       " characters but only ", len(inputs%imp_positions_str), &
                                       " are kept."
                print *, "  The list would be truncated in silence and the run would use fewer"
                print *, "  impurities than requested. Shorten the list."
                ierr = ERROR_INVALID_INPUT
                return
            end if
        else
            call print_missing_group_note('potential', filename)
        end if

        if (namelist_group_present(records, 'scf')) then
            io_msg = ''
            read(records, nml=scf, iostat=io_stat, iomsg=io_msg)
            call check_namelist_status('scf', filename, io_stat, io_msg, ierr)
            if (ierr /= ERROR_SUCCESS) return
        else
            call print_missing_group_note('scf', filename)
        end if

        if (namelist_group_present(records, 'output')) then
            io_msg = ''
            read(records, nml=output, iostat=io_stat, iomsg=io_msg)
            call check_namelist_status('output', filename, io_stat, io_msg, ierr)
            if (ierr /= ERROR_SUCCESS) return
        else
            call print_missing_group_note('output', filename)
        end if

        ! Update inputs structure
        inputs%L = L
        inputs%Nup = Nup
        inputs%Ndown = Ndown
        inputs%U = U
        inputs%bc_type = bc
        inputs%phase = phase
        inputs%table_dir = table_dir

        inputs%potential_type = potential_type
        inputs%V0 = V0
        inputs%pot_center = pot_center
        inputs%pot_width = pot_width
        inputs%concentration = concentration
        inputs%pot_seed = pot_seed
        ! The explicit substring is the whole destination field: the guard after
        ! the &potential read has already rejected anything longer, so this copy
        ! cannot lose characters. Written this way so that the (correct)
        ! -Wcharacter-truncation warning for a 4000 -> 500 assignment does not
        ! have to be silenced globally.
        inputs%imp_positions_str = imp_positions_str(1:len(inputs%imp_positions_str))
        inputs%disorder_strength = disorder_strength
        inputs%spring_constant = spring_constant
        inputs%aah_lambda = aah_lambda
        inputs%aah_beta = aah_beta
        inputs%aah_phi = aah_phi
        inputs%position = position
        inputs%width = width
        inputs%barrier_width = barrier_width
        inputs%well_depth = well_depth
        inputs%well_width = well_width
        inputs%position1 = position1
        inputs%width1 = width1
        inputs%position2 = position2
        inputs%width2 = width2
        
        inputs%max_iter = max_iter
        inputs%density_tol = density_tol
        inputs%energy_tol = energy_tol
        inputs%potential_tol = potential_tol
        inputs%mixing_alpha = mixing_alpha
        inputs%verbose = verbose
        inputs%store_history = store_history
        inputs%use_adaptive_mixing = use_adaptive_mixing
        inputs%xc_smoothing_width = xc_smoothing_width
        
        inputs%output_prefix = output_prefix
        inputs%save_density = save_density
        inputs%save_eigenvalues = save_eigenvalues
        inputs%save_wavefunction = save_wavefunction
        
    end subroutine read_namelist_file

    !> Slurp a text file into an array of complete, fixed-length records
    !!
    !! The array is meant to be used as an INTERNAL file for the namelist reads
    !! (see read_namelist_file): unlike an external unit, an internal file has no
    !! notion of a "newline", so the last group of a file whose last line is not
    !! newline-terminated is read exactly like any other.
    !!
    !! Records are read with a non-advancing loop and concatenated, so a line of
    !! any length is preserved in full: truncating a line would silently change
    !! the meaning of a long value (e.g. imp_positions_str) instead of failing.
    !! The length of the returned array is the length of the longest line, and
    !! shorter lines are blank-padded, which is what a formatted read expects.
    !!
    !! An empty file yields a single blank record, so that the namelist reads
    !! report "group absent" (the normal, non-fatal path) rather than failing on
    !! a zero-sized internal file.
    !!
    !! @param[in]  filename File to read
    !! @param[out] records  One element per line of the file, blank-padded
    !! @param[out] ierr     ERROR_SUCCESS, ERROR_FILE_NOT_FOUND or
    !!                      ERROR_INVALID_INPUT (unreadable content)
    subroutine read_input_records(filename, records, ierr)
        character(len=*), intent(in) :: filename
        character(len=:), allocatable, intent(out) :: records(:)
        integer, intent(out) :: ierr

        integer, parameter :: CHUNK_LEN = 512
        character(len=CHUNK_LEN) :: chunk
        integer :: io_unit, io_stat, n_read, n_rec, max_len, cur_len, k

        ierr = ERROR_SUCCESS

        open(newunit=io_unit, file=filename, status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            print *, "ERROR: Could not open file: ", trim(filename)
            ierr = ERROR_FILE_NOT_FOUND
            return
        end if

        ! First pass: count the records and measure the longest one.
        n_rec = 0
        max_len = 0
        cur_len = 0
        do
            read(io_unit, '(A)', advance='no', size=n_read, iostat=io_stat) chunk
            cur_len = cur_len + max(n_read, 0)

            if (io_stat == IOSTAT_EOR) then
                n_rec = n_rec + 1
                max_len = max(max_len, cur_len)
                cur_len = 0
            else if (io_stat == IOSTAT_END) then
                ! A last line with no terminating newline is reported as a
                ! complete record by gfortran, but the standard allows the
                ! end-of-file condition here, so the pending characters are
                ! counted as one more record.
                if (cur_len > 0) then
                    n_rec = n_rec + 1
                    max_len = max(max_len, cur_len)
                end if
                exit
            else if (io_stat /= 0) then
                print *, "ERROR: Could not read file: ", trim(filename)
                close(io_unit)
                ierr = ERROR_INVALID_INPUT
                return
            end if
            ! io_stat == 0: the record is longer than CHUNK_LEN, keep going.
        end do

        rewind(io_unit)

        allocate(character(len=max(max_len, 1)) :: records(max(n_rec, 1)))
        records = ''

        ! Second pass: each record fits in one element, so a plain advancing
        ! read transfers it whole.  A failed read here means the contents no
        ! longer match the record count established in the first pass; leaving
        ! the remaining entries blank would make later namelist groups look
        ! absent and silently retain their defaults.
        do k = 1, n_rec
            read(io_unit, '(A)', iostat=io_stat) records(k)
            if (io_stat /= 0) then
                print *, "ERROR: Could not read file: ", trim(filename)
                close(io_unit)
                ierr = ERROR_INVALID_INPUT
                return
            end if
        end do

        close(io_unit)
    end subroutine read_input_records

    !> Decide whether a namelist `read` status is acceptable, and report it
    !!
    !! Called ONLY for a group whose header was found in the file, so every
    !! non-zero status here means "the group is there and it is broken". This
    !! is what lets the routine be strict: the ambiguity that used to force a
    !! second interpretation of `iostat` lives in the compiler's reporting of an
    !! ABSENT group, and read_namelist_file no longer relies on it.
    !!
    !! Measured behaviour of GNU Fortran 16.1 on a namelist `read`:
    !!
    !! | situation                                | iostat | iomsg                          |
    !! |------------------------------------------|--------|--------------------------------|
    !! | group read normally                      |      0 | (undefined)                    |
    !! | group present but empty (`&g /`)         |      0 | (undefined)                    |
    !! | unknown/misspelled key inside the group  |   5010 | "Cannot match namelist object  |
    !! |                                          |        |  name <key>"                   |
    !! | malformed value (`L = abc`)              |   5010 | "Cannot match namelist object  |
    !! |                                          |        |  name abc"                     |
    !! | malformed value (`U = 1.2.3`), external  |     -1 | "End of file"                  |
    !! | malformed value (`U = 1.2.3`), internal  |   5010 | "Cannot match namelist object  |
    !! |                                          |        |  name .2.3"                    |
    !! | group not terminated by `/`              |     -1 | "End of file"                  |
    !! | group absent, external unit              |     -1 | "End of file"                  |
    !! | group absent, internal file              |      0 | (undefined)                    |
    !! | last line not newline-terminated,        |        |                                |
    !! | external unit only                       |     -1 | "End of file"                  |
    !!
    !! The last row is why read_namelist_file reads from an internal file: on an
    !! external unit a missing final newline is indistinguishable from a broken
    !! group, and a perfectly valid input file was rejected because of it.
    !!
    !! @param[in]  group    Name of the namelist group, without the leading '&'
    !! @param[in]  filename Input file being read (for the messages only)
    !! @param[in]  io_stat  iostat returned by the namelist read
    !! @param[in]  io_msg   iomsg returned by the namelist read
    !! @param[out] ierr     ERROR_SUCCESS or ERROR_INVALID_INPUT
    subroutine check_namelist_status(group, filename, io_stat, io_msg, ierr)
        character(len=*), intent(in) :: group
        character(len=*), intent(in) :: filename
        integer, intent(in) :: io_stat
        character(len=*), intent(in) :: io_msg
        integer, intent(out) :: ierr

        ierr = ERROR_SUCCESS
        if (io_stat == 0) return

        print *, "ERROR reading &" // group // " in " // trim(filename) // ": " // trim(io_msg)
        if (io_stat < 0) then
            print *, "  The group is in the file but the read stopped before its end."
            print *, "  Typical causes: a malformed value (e.g. U = 1.2.3, L = abc),"
            print *, "  an unquoted string, or a group with no closing '/'."
        end if
        call print_namelist_hint(group, io_msg)
        ierr = ERROR_INVALID_INPUT
    end subroutine check_namelist_status

    !> Report a namelist group that is not in the input file at all
    !!
    !! Not an error: the four groups are optional by design (examples/
    !! input_minimal.txt has &system only) and an absent group leaves every one
    !! of its keys at the default. It is printed all the same, because a group
    !! that the user believes is being read - a typo in the header, a group
    !! commented out - looks exactly like this.
    !!
    !! @param[in] group    Name of the namelist group, without the leading '&'
    !! @param[in] filename Input file being read
    subroutine print_missing_group_note(group, filename)
        character(len=*), intent(in) :: group
        character(len=*), intent(in) :: filename

        print *, "NOTE: group &" // group // " not found in " // trim(filename) // &
                 "; its parameters keep their default values."
    end subroutine print_missing_group_note

    !> Print a hint for a failed namelist read, naming obsolete keys explicitly
    !!
    !! An input file written for an older version of the code fails on a key
    !! that was deliberately removed, and the compiler message ("Cannot match
    !! namelist object name distribution") gives the user no way to tell a typo
    !! from a removal. The obsolete keys are therefore listed by name, with the
    !! replacement, so the file can be migrated instead of guessed at.
    !!
    !! @param[in] group  Namelist group that failed
    !! @param[in] io_msg iomsg from the failed read (searched for obsolete keys)
    subroutine print_namelist_hint(group, io_msg)
        character(len=*), intent(in) :: group
        character(len=*), intent(in) :: io_msg

        character(len=len(io_msg)) :: msg_lower

        msg_lower = to_lower(io_msg)

        if (index(msg_lower, 'distribution') > 0) then
            print *, "  NOTE: the key 'distribution' was REMOVED from &potential."
            print *, "        The disorder distribution is selected by potential_type alone:"
            print *, "        potential_type = 'random_uniform' or 'random_gaussian'."
            print *, "        Delete the 'distribution' line from the input file."
            return
        end if

        print *, "  Check &" // group // " for a misspelled key, for a key that belongs"
        print *, "  to another group, or for a key removed in a newer version of the code."
    end subroutine print_namelist_hint

    !> True if the lines of an input file contain the header of a namelist group
    !!
    !! Used to decide, before the read, whether a group is in the file at all:
    !! the compiler cannot be asked, since it reports an absent group with the
    !! same status as a well-formed one (internal file) or as a broken one
    !! (external unit) - see check_namelist_status. A header is a '&' or a '$'
    !! that is not inside a quoted string and not in a comment, followed by the
    !! group name (namelist group names are case-insensitive) and then a blank,
    !! a comma, a slash, an '=' or the end of the line.
    !!
    !! BOTH header characters must be accepted: gfortran also reads the legacy
    !! form `$system ... $end` (measured: iostat = 0 and every value assigned).
    !! Recognising only '&' made this function return .false. for a file the
    !! reader handles perfectly, and since the decision now gates the read that
    !! turned the whole group into a silent fallback to the defaults - the exact
    !! failure mode the iostat checking exists to eliminate, and in its worst
    !! direction, because the run then completes with plausible numbers for a
    !! system nobody asked for. `$end` is the terminator of that form, so the
    !! name `end` is never reported as a group.
    !!
    !! The scan works on the records produced by read_input_records instead of
    !! re-opening the file: those are the exact lines the namelist read sees
    !! (including a last line with no terminating newline), and it avoids a
    !! second, racy pass over the file.
    !!
    !! The decision now gates the read, so BOTH kinds of mistake are costly and
    !! both are guarded against:
    !!
    !! * a false positive would make an absent group abort the run: hence '&'
    !!   inside a quoted value (`output_prefix = '&potential'`) and everything
    !!   after a '!' comment marker are skipped;
    !! * a false negative would make a group that IS in the file be ignored in
    !!   silence, with every one of its keys quietly back at the default: hence
    !!   both header characters are accepted, tabs and carriage returns are
    !!   treated as blanks (`adjustl` does not move them, and a record coming
    !!   from a CRLF file may still carry the CR) and the whole line is scanned,
    !!   not just its first character, so a header that does not start the line
    !!   - `&system L=8 / &scf max_iter=5 /` - is still found.
    !!
    !! @param[in] records Lines of the input file, blank-padded
    !! @param[in] group   Group name, without the leading '&'
    !! @return    .true. if a header for `group` was found
    pure function namelist_group_present(records, group) result(is_present)
        character(len=*), intent(in) :: records(:)
        character(len=*), intent(in) :: group
        logical :: is_present

        character, parameter :: TAB = char(9)
        character, parameter :: CR = char(13)
        integer, parameter :: MAX_NAME = 64
        character(len=len(records)) :: line
        character(len=MAX_NAME) :: name
        character :: ch, quote, header
        integer :: k, i, j, n, n_line

        is_present = .false.

        do k = 1, size(records)
            line = records(k)
            do i = 1, len(line)
                if (line(i:i) == TAB .or. line(i:i) == CR) line(i:i) = ' '
            end do

            n_line = len_trim(line)
            i = 1
            do while (i <= n_line)
                ch = line(i:i)

                if (ch == '!') exit  ! comment: the rest of the line is prose

                if (ch == '''' .or. ch == '"') then
                    ! Skip the quoted string, closing quote included.
                    quote = ch
                    i = i + 1
                    do while (i <= n_line)
                        if (line(i:i) == quote) exit
                        i = i + 1
                    end do
                    i = i + 1
                    cycle
                end if

                if (ch /= '&' .and. ch /= '$') then
                    i = i + 1
                    cycle
                end if

                ! Candidate header: collect the group name that follows.
                header = ch
                n = 0
                j = i + 1
                do while (j <= n_line)
                    ch = line(j:j)
                    if (ch == ' ' .or. ch == ',' .or. ch == '/' .or. &
                        ch == '=' .or. ch == '!') exit
                    if (n >= MAX_NAME) exit
                    n = n + 1
                    name(n:n) = ch
                    j = j + 1
                end do

                if (n > 0) then
                    ! `$end` closes a group of the legacy form; it is a
                    ! terminator, never the name of a group.
                    if (.not. (header == '$' .and. to_lower(name(1:n)) == 'end')) then
                        if (to_lower(name(1:n)) == to_lower(group)) then
                            is_present = .true.
                            return
                        end if
                    end if
                end if

                i = max(j, i + 1)
            end do
        end do
    end function namelist_group_present

    !> Lowercase copy of a string (ASCII only)
    !!
    !! @param[in] str String to fold
    !! @return    Same string with 'A'-'Z' replaced by 'a'-'z'
    pure function to_lower(str) result(folded)
        character(len=*), intent(in) :: str
        character(len=len(str)) :: folded
        integer :: i, code

        do i = 1, len(str)
            code = iachar(str(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                folded(i:i) = achar(code - iachar('A') + iachar('a'))
            else
                folded(i:i) = str(i:i)
            end if
        end do
    end function to_lower

    !> Validate input parameters
    !!
    !! @param[in]  inputs Input parameters
    !! @param[out] ierr   Error code (0 = success)
    subroutine validate_inputs(inputs, ierr)
        type(input_params_t), intent(in) :: inputs
        integer, intent(out) :: ierr
        
        ierr = ERROR_SUCCESS
        
        ! Validate system parameters
        if (inputs%L <= 0) then
            print *, "ERROR: L must be positive, got L =", inputs%L
            ierr = ERROR_INVALID_INPUT
            return
        end if
        
        if (inputs%Nup < 0 .or. inputs%Ndown < 0) then
            print *, "ERROR: Nup and Ndown must be non-negative"
            ierr = ERROR_INVALID_INPUT
            return
        end if
        
        ! Pauli exclusion per spin channel: each spin channel has exactly L
        ! orbitals, so 0 <= Nup <= L and 0 <= Ndown <= L. The total bound
        ! Nup + Ndown <= 2L is a consequence of these two.
        if (inputs%Nup > inputs%L) then
            print *, "ERROR: Nup cannot exceed L (only L spin-up orbitals available)"
            print *, "  Nup =", inputs%Nup, ", L =", inputs%L
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (inputs%Ndown > inputs%L) then
            print *, "ERROR: Ndown cannot exceed L (only L spin-down orbitals available)"
            print *, "  Ndown =", inputs%Ndown, ", L =", inputs%L
            ierr = ERROR_INVALID_INPUT
            return
        end if
        
        ! Validate boundary conditions
        select case (trim(inputs%bc_type))
        case ('open', 'periodic', 'twisted')
            ! Valid
        case default
            print *, "ERROR: Invalid boundary condition: ", trim(inputs%bc_type)
            print *, "  Valid options: open, periodic, twisted"
            ierr = ERROR_INVALID_INPUT
            return
        end select

        ! The twist angle is given in units of π, so the full physical range of
        ! the Aharonov-Bohm phase is [0, 2). After the conversion to radians in
        ! main.f90 this is exactly the interval [0, 2π) demanded by
        ! validate_bc_parameters; rejecting it here turns a cryptic
        ! ERROR_OUT_OF_BOUNDS raised deep inside the Hamiltonian builder into a
        ! message that names the offending input.
        !
        ! NOTE: this is a DELIBERATE divergence from the C++ reference, which
        ! validates only `bc` and would silently accept a negative phase or a
        ! phase above 2π (`original/lsda_interface.cc`).
        if (trim(inputs%bc_type) == 'twisted') then
            if (.not. ieee_is_finite(inputs%phase) .or. &
                inputs%phase < 0.0_dp .or. inputs%phase >= 2.0_dp) then
                print *, "ERROR: phase must be in [0, 2) in units of pi (i.e. [0, 2*pi) radians), got:", &
                         inputs%phase
                ierr = ERROR_INVALID_INPUT
                return
            end if
        end if

        ! Validate concentration for impurity potential
        if (trim(inputs%potential_type) == 'impurity') then
            if (inputs%concentration <= 0.0_dp .or. &
                inputs%concentration > 100.0_dp) then
                print *, "ERROR: concentration must be in (0, 100], got:", &
                                                        inputs%concentration
                ierr = ERROR_INVALID_INPUT
                return
            end if
        end if
        
        ! Validate SCF parameters
        if (inputs%max_iter <= 0) then
            print *, "ERROR: max_iter must be positive"
            ierr = ERROR_INVALID_INPUT
            return
        end if
        
        ! potential_tol and energy_tol are BOTH convergence criteria (the SCF
        ! stops on residual_V < potential_tol AND |dE| < energy_tol*max(1,|E|)).
        ! A null, negative or NaN tolerance can never be met, so the run would be
        ! condemned to exhaust max_iter and report a convergence failure instead
        ! of an invalid input. An infinite tolerance is worse: its comparison is
        ! always true for finite residuals/energies, which silently switches that
        ! criterion off and lets the cycle declare convergence on the other one.
        if (.not. ieee_is_finite(inputs%potential_tol) .or. &
            inputs%potential_tol <= 0.0_dp) then
            print *, "ERROR: potential_tol must be finite and positive, got:", inputs%potential_tol
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (.not. ieee_is_finite(inputs%energy_tol) .or. &
            inputs%energy_tol <= 0.0_dp) then
            print *, "ERROR: energy_tol must be finite and positive, got:", inputs%energy_tol
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (inputs%mixing_alpha <= 0.0_dp .or. inputs%mixing_alpha > 1.0_dp) then
            print *, "ERROR: mixing_alpha must be in (0, 1]"
            ierr = ERROR_INVALID_INPUT
            return
        end if

        ! Smoothing of the V_xc discontinuity at n = 1: 0 disables it (exact C++
        ! parity); the window [1-w, 1+w] must stay inside the physical range of n.
        ! NaN would pass both range comparisons and then fail SILENTLY: get_vxc
        ! evaluates `w > 0` as false and quietly takes the unsmoothed branch,
        ! so the user asks for smoothing and gets none. Both infinities are
        ! already rejected by the range checks (+Inf >= 1, -Inf < 0).
        if (ieee_is_nan(inputs%xc_smoothing_width) .or. &
            inputs%xc_smoothing_width < 0.0_dp .or. &
            inputs%xc_smoothing_width >= XC_SMOOTHING_WIDTH_MAX) then
            print *, "ERROR: xc_smoothing_width must be in [0, 1), got:", inputs%xc_smoothing_width
            ierr = ERROR_INVALID_INPUT
            return
        end if
        
    end subroutine validate_inputs

    !> Convert input_params_t to system_params_t
    !!
    !! This is the SINGLE place where the twist angle changes unit:
    !! `inputs%phase` is in units of π (the user-facing convention, identical to
    !! the C++ reference, whose prompt says "in units of Pi" and which applies
    !! `phase *= M_PI` right after reading it) and `sys_params%phase` is in
    !! radians, the unit expected by `apply_boundary_conditions_complex` and
    !! validated as [0, 2π) by `validate_bc_parameters`.
    !!
    !! The multiplication used to live in the caller (`app/main.f90`), which
    !! left this routine returning a `system_params_t` whose `phase` was in the
    !! wrong unit for every one of its consumers. With a single caller the
    !! contract closed by accident; any new executable, example or test that
    !! converted the inputs and handed the result to the Kohn-Sham cycle would
    !! have run with θ = phase rad instead of θ = phase·π, and nothing
    !! downstream can tell 0.5 from 1.5708.
    !!
    !! @param[in]  inputs     Input parameters (phase in units of π)
    !! @param[out] sys_params System parameters (phase in radians)
    !! @param[out] ierr       Error code
    subroutine convert_to_system_params(inputs, sys_params, ierr)
        type(input_params_t), intent(in) :: inputs
        type(system_params_t), intent(out) :: sys_params
        integer, intent(out) :: ierr
        
        ierr = ERROR_SUCCESS
        
        sys_params%L = inputs%L
        sys_params%Nup = inputs%Nup
        sys_params%Ndown = inputs%Ndown
        sys_params%U = inputs%U
        sys_params%phase = inputs%phase * PI
        
        ! Convert BC string to integer
        select case (trim(inputs%bc_type))
        case ('open')
            sys_params%bc = BC_OPEN
        case ('periodic')
            sys_params%bc = BC_PERIODIC
        case ('twisted')
            sys_params%bc = BC_TWISTED
        case default
            ierr = ERROR_INVALID_INPUT
        end select
        
    end subroutine convert_to_system_params

    !> Convert input_params_t to scf_params_t
    !!
    !! @param[in]  inputs     Input parameters
    !! @param[out] scf_params SCF parameters
    subroutine convert_to_scf_params(inputs, scf_params)
        type(input_params_t), intent(in) :: inputs
        type(scf_params_t), intent(out) :: scf_params
        
        scf_params%max_iter = inputs%max_iter
        scf_params%density_tol = inputs%density_tol
        scf_params%energy_tol = inputs%energy_tol
        scf_params%potential_tol = inputs%potential_tol
        scf_params%mixing_alpha = inputs%mixing_alpha
        scf_params%verbose = inputs%verbose
        scf_params%store_history = inputs%store_history
        scf_params%use_adaptive_mixing = inputs%use_adaptive_mixing
        
    end subroutine convert_to_scf_params

    !> Lattice bounds of a single rectangular barrier of exactly `width` sites
    !!
    !! The barrier is centred on `position` and spans `[i_start, i_end]` with
    !! `i_end - i_start + 1 == width` for **every** width, even or odd.
    !!
    !! The previous expression `i_start = position - width/2`,
    !! `i_end = position + width/2` relied on integer division on both sides and
    !! therefore produced `2*(width/2) + 1` sites: correct for odd widths, but
    !! `width + 1` sites for even widths (width = 4 gave 5 occupied sites). Since
    !! the barrier width controls the tunnelling probability exponentially, an
    !! off-by-one site is a physics error, not a cosmetic one.
    !!
    !! For even widths the extra site is taken on the left, i.e. the barrier
    !! covers `[position - width/2, position + width/2 - 1]`, which keeps
    !! `position` inside the barrier.
    !!
    !! @param[in]  position Centre site of the barrier (1-indexed)
    !! @param[in]  width    Number of sites covered by the barrier (must be > 0)
    !! @param[out] i_start  First site of the barrier
    !! @param[out] i_end    Last site of the barrier
    subroutine barrier_single_bounds(position, width, i_start, i_end)
        integer, intent(in) :: position, width
        integer, intent(out) :: i_start, i_end

        i_start = position - width / 2
        i_end = i_start + width - 1
    end subroutine barrier_single_bounds

    !> Parse a comma/blank separated list of lattice positions into integers
    !!
    !! Used for the `imp_positions_str` field of the `&potential` namelist, e.g.
    !! `imp_positions_str = '10, 25, 40, 55'`.
    !!
    !! The parser is deliberately strict: a malformed list must produce a clear
    !! error instead of a silently truncated list of impurities, because a
    !! truncated list still yields a perfectly plausible - and wrong -
    !! calculation. The following are all rejected with ERROR_INVALID_INPUT:
    !!
    !! - an empty (or blank-only) string;
    !! - a leading or trailing comma, or two consecutive commas (empty field);
    !! - a token that is not a plain optionally signed integer (`'12a'`, `'1.5'`);
    !! - a position outside `[1, max_value]`;
    !! - a repeated position.
    !!
    !! Blanks around the commas are irrelevant, and a run of blanks without any
    !! comma also separates two tokens (`'10 25'` is the same as `'10, 25'`).
    !!
    !! @param[in]  str        String to parse
    !! @param[in]  max_value  Upper bound for the accepted positions (normally L)
    !! @param[out] values     Parsed positions, in the order given (size = count).
    !!                        Left unallocated when ierr /= ERROR_SUCCESS.
    !! @param[out] ierr       ERROR_SUCCESS or ERROR_INVALID_INPUT
    subroutine parse_int_list(str, max_value, values, ierr)
        character(len=*), intent(in) :: str
        integer, intent(in) :: max_value
        integer, allocatable, intent(out) :: values(:)
        integer, intent(out) :: ierr

        character, parameter :: TAB = char(9)
        integer, allocatable :: work(:)
        integer :: n, i, j, n_tok, n_comma, value, io_stat
        character(len=len(str)) :: token
        character :: ch

        ierr = ERROR_SUCCESS

        n = len_trim(str)
        if (n == 0) then
            print *, "ERROR: imp_positions_str is empty; give at least one site, e.g. '10, 25, 40'"
            ierr = ERROR_INVALID_INPUT
            return
        end if

        ! Upper bound on the number of tokens: one token needs at least one
        ! character, so n tokens is the worst case.
        allocate(work(n))
        n_tok = 0
        n_comma = 0
        i = 1

        do while (i <= n)
            ch = str(i:i)

            if (ch == ',') then
                n_comma = n_comma + 1
                if (n_comma > 1) then
                    print *, "ERROR: empty field in imp_positions_str (two consecutive separators): ", &
                             trim(str)
                    ierr = ERROR_INVALID_INPUT
                    return
                end if
                i = i + 1
                cycle
            end if

            if (ch == ' ' .or. ch == TAB) then
                i = i + 1
                cycle
            end if

            ! Start of a token
            if (n_tok == 0 .and. n_comma > 0) then
                print *, "ERROR: imp_positions_str starts with a separator: ", trim(str)
                ierr = ERROR_INVALID_INPUT
                return
            end if

            j = i
            do while (j <= n)
                ch = str(j:j)
                if (ch == ',' .or. ch == ' ' .or. ch == TAB) exit
                j = j + 1
            end do

            token = str(i:j-1)
            if (.not. is_plain_integer(token(1:j-i))) then
                print *, "ERROR: '", token(1:j-i), "' is not an integer in imp_positions_str: ", trim(str)
                ierr = ERROR_INVALID_INPUT
                return
            end if

            read(token(1:j-i), *, iostat=io_stat) value
            if (io_stat /= 0) then
                print *, "ERROR: could not read the integer '", token(1:j-i), "' in imp_positions_str"
                ierr = ERROR_INVALID_INPUT
                return
            end if

            if (value < 1 .or. value > max_value) then
                print *, "ERROR: impurity site out of range [1, ", max_value, "] in imp_positions_str:", value
                ierr = ERROR_INVALID_INPUT
                return
            end if

            if (n_tok > 0) then
                if (any(work(1:n_tok) == value)) then
                    print *, "ERROR: repeated impurity site in imp_positions_str:", value
                    ierr = ERROR_INVALID_INPUT
                    return
                end if
            end if

            n_tok = n_tok + 1
            work(n_tok) = value
            n_comma = 0
            i = j
        end do

        if (n_comma > 0) then
            print *, "ERROR: imp_positions_str ends with a separator: ", trim(str)
            ierr = ERROR_INVALID_INPUT
            return
        end if

        ! Unreachable for a non-blank string (the blank-only case is caught by
        ! len_trim above), kept as a guard against future changes of the loop.
        if (n_tok == 0) then
            print *, "ERROR: no impurity site found in imp_positions_str: ", trim(str)
            ierr = ERROR_INVALID_INPUT
            return
        end if

        allocate(values(n_tok))
        values = work(1:n_tok)
    end subroutine parse_int_list

    !> True if the token is an optionally signed sequence of decimal digits
    !!
    !! Needed because list-directed `read` happily accepts `'1.5'` (as 1) and
    !! stops at the first offending character of `'12a'`, which would let a
    !! typo through as a valid - but different - lattice site.
    !!
    !! @param[in] token Token to inspect (must not contain blanks)
    !! @return    .true. if the token is a plain integer literal
    pure function is_plain_integer(token) result(ok)
        character(len=*), intent(in) :: token
        logical :: ok
        integer :: k, first

        ok = .false.
        if (len(token) == 0) return

        first = 1
        if (token(1:1) == '+' .or. token(1:1) == '-') first = 2
        if (first > len(token)) return

        do k = first, len(token)
            if (token(k:k) < '0' .or. token(k:k) > '9') return
        end do

        ok = .true.
    end function is_plain_integer

    !> Parse integer argument from command line
    subroutine parse_integer_arg(current_idx, nargs, value, ierr)
        integer, intent(in) :: current_idx, nargs
        integer, intent(out) :: value
        integer, intent(out) :: ierr

        character(len=256) :: arg, arg_name
        integer :: io_stat

        if (current_idx + 1 > nargs) then
            call get_command_argument(current_idx, arg_name)
            print *, "ERROR: ", trim(arg_name), " requires a value"
            ierr = ERROR_INVALID_INPUT
            return
        end if

        call get_command_argument(current_idx + 1, arg)
        read(arg, *, iostat=io_stat) value

        if (io_stat /= 0) then
            call get_command_argument(current_idx, arg_name)
            print *, "ERROR: Invalid integer for ", trim(arg_name), ": ", trim(arg)
            ierr = ERROR_INVALID_INPUT
        else
            ierr = ERROR_SUCCESS
        end if
    end subroutine parse_integer_arg

    !> Parse real argument from command line
    subroutine parse_real_arg(current_idx, nargs, value, ierr)
        integer, intent(in) :: current_idx, nargs
        real(dp), intent(out) :: value
        integer, intent(out) :: ierr

        character(len=256) :: arg, arg_name
        integer :: io_stat

        if (current_idx + 1 > nargs) then
            call get_command_argument(current_idx, arg_name)
            print *, "ERROR: ", trim(arg_name), " requires a value"
            ierr = ERROR_INVALID_INPUT
            return
        end if

        call get_command_argument(current_idx + 1, arg)
        read(arg, *, iostat=io_stat) value

        if (io_stat /= 0) then
            call get_command_argument(current_idx, arg_name)
            print *, "ERROR: Invalid real number for ", trim(arg_name), ": ", trim(arg)
            ierr = ERROR_INVALID_INPUT
        else
            ierr = ERROR_SUCCESS
        end if
    end subroutine parse_real_arg

    !> Print help message
    subroutine print_help()
        print '(A)', ""
        print '(A)', "LSDA-Hubbard: 1D Hubbard Model DFT-LSDA Calculator"
        print '(A)', "=================================================="
        print '(A)', ""
        print '(A)', "Usage:"
        print '(A)', "  fpm run lsdaks -- --input <file>     # Read from namelist file"
        print '(A)', "  fpm run lsdaks -- [options]          # Command line mode"
        print '(A)', ""
        print '(A)', "System Parameters:"
        print '(A)', "  --L <int>         Number of lattice sites (default: 10)"
        print '(A)', "  --Nup <int>       Number of spin-up electrons (default: 5)"
        print '(A)', "  --Ndown <int>     Number of spin-down electrons (default: 5)"
        print '(A)', "  --U <real>        Hubbard interaction strength (default: 4.0)"
        print '(A)', "  --bc <type>       Boundary conditions: open/periodic/twisted (default: periodic)"
        print '(A)', "  --phase <real>    Twist angle for twisted BC, in units of pi, in [0, 2) (default: 0.0)"
        print '(A)', ""
        print '(A)', "Potential:"
        print '(A)', "  --potential <type>    Type: uniform/harmonic/impurity/... (default: uniform)"
        print '(A)', "  --V0 <real>           Potential strength (default: 0.0)"
        print '(A)', "  --concentration <real> Impurity concentration in % (default: 50.0)"
        print '(A)', ""
        print '(A)', "Output:"
        print '(A)', "  --verbose         Print SCF progress (default: on)"
        print '(A)', "  --quiet           Suppress SCF progress"
        print '(A)', ""
        print '(A)', "Other:"
        print '(A)', "  --help, -h        Show this help message"
        print '(A)', ""
        print '(A)', "Examples:"
        print '(A)', "  fpm run lsdaks -- --L 10 --Nup 5 --Ndown 5 --U 4.0"
        print '(A)', "  fpm run lsdaks -- --input my_simulation.txt"
        print '(A)', ""
    end subroutine print_help

end module input_parser
