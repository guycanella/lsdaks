!> Module for parsing input parameters from command line and namelist files
!!
!! Supports two input modes:
!! 1. Namelist file: fpm run lsdaks -- --input input.txt
!! 2. Command line:  fpm run lsdaks -- --L 10 --Nup 5 --Ndown 5 --U 4.0
!!
!! Priority: Command line arguments override namelist values
module input_parser
    use lsda_constants, only: dp, ITER_MAX, SCF_DENSITY_TOL, &
                                    SCF_ENERGY_TOL, MIX_ALPHA
    use lsda_types, only: system_params_t
    use kohn_sham_cycle, only: scf_params_t
    use boundary_conditions, only: BC_OPEN, BC_PERIODIC, BC_TWISTED
    use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT, &
                                            ERROR_FILE_NOT_FOUND
    implicit none
    private

    !> Input parameters type
    type, public :: input_params_t
        ! System parameters
        integer :: L = 10
        integer :: Nup = 5
        integer :: Ndown = 5
        real(dp) :: U = 4.0_dp
        character(len=20) :: bc_type = 'periodic'
        real(dp) :: phase = 0.0_dp
        
        ! Potential parameters
        character(len=20) :: potential_type = 'uniform'
        real(dp) :: V0 = 0.0_dp
        real(dp) :: pot_center = 0.0_dp
        real(dp) :: pot_width = 1.0_dp
        
        ! SCF parameters
        integer :: max_iter = ITER_MAX
        real(dp) :: density_tol = SCF_DENSITY_TOL
        real(dp) :: energy_tol = SCF_ENERGY_TOL
        real(dp) :: mixing_alpha = MIX_ALPHA
        logical :: verbose = .true.
        logical :: store_history = .true.
        
        ! Output parameters
        character(len=100) :: output_prefix = 'lsda_output'
        logical :: save_density = .true.
        logical :: save_eigenvalues = .true.
        logical :: save_wavefunction = .false.
    end type input_params_t

    public :: parse_inputs
    public :: validate_inputs
    public :: convert_to_system_params
    public :: convert_to_scf_params

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
        
        inputs = input_params_t()
        
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
                    
                case ('--verbose')
                    inputs%verbose = .true.
                    i = i + 1
                    
                case ('--quiet')
                    inputs%verbose = .false.
                    i = i + 1
                    
                case ('--help', '-h')
                    call print_help()
                    ierr = -1  ! Special code: help requested, exit gracefully
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
    !! @param[in]    filename Input file path
    !! @param[inout] inputs   Input parameters (updated from file)
    !! @param[out]   ierr     Error code
    subroutine read_namelist_file(filename, inputs, ierr)
        character(len=*), intent(in) :: filename
        type(input_params_t), intent(inout) :: inputs
        integer, intent(out) :: ierr
        
        integer :: io_unit, io_stat
        logical :: file_exists
        
        ! System namelist variables
        integer :: L, Nup, Ndown
        real(dp) :: U, phase
        character(len=20) :: bc
        
        ! Potential namelist variables
        character(len=20) :: potential_type
        real(dp) :: V0, pot_center, pot_width
        
        ! SCF namelist variables
        integer :: max_iter
        real(dp) :: density_tol, energy_tol, mixing_alpha
        logical :: verbose, store_history
        
        ! Output namelist variables
        character(len=100) :: output_prefix
        logical :: save_density, save_eigenvalues, save_wavefunction
        
        namelist /system/ L, Nup, Ndown, U, bc, phase
        namelist /potential/ potential_type, V0, pot_center, pot_width
        namelist /scf/ max_iter, density_tol, energy_tol, mixing_alpha, verbose, store_history
        namelist /output/ output_prefix, save_density, save_eigenvalues, save_wavefunction
        
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
        
        potential_type = inputs%potential_type
        V0 = inputs%V0
        pot_center = inputs%pot_center
        pot_width = inputs%pot_width
        
        max_iter = inputs%max_iter
        density_tol = inputs%density_tol
        energy_tol = inputs%energy_tol
        mixing_alpha = inputs%mixing_alpha
        verbose = inputs%verbose
        store_history = inputs%store_history
        
        output_prefix = inputs%output_prefix
        save_density = inputs%save_density
        save_eigenvalues = inputs%save_eigenvalues
        save_wavefunction = inputs%save_wavefunction
        
        ! Open and read file
        open(newunit=io_unit, file=filename, status='old', iostat=io_stat)
        if (io_stat /= 0) then
            print *, "ERROR: Could not open file: ", trim(filename)
            ierr = ERROR_FILE_NOT_FOUND
            return
        end if
        
        ! Read namelists (silently skip missing sections)
        read(io_unit, nml=system, iostat=io_stat)
        rewind(io_unit)
        
        read(io_unit, nml=potential, iostat=io_stat)
        rewind(io_unit)
        
        read(io_unit, nml=scf, iostat=io_stat)
        rewind(io_unit)
        
        read(io_unit, nml=output, iostat=io_stat)
        
        close(io_unit)
        
        ! Update inputs structure
        inputs%L = L
        inputs%Nup = Nup
        inputs%Ndown = Ndown
        inputs%U = U
        inputs%bc_type = bc
        inputs%phase = phase
        
        inputs%potential_type = potential_type
        inputs%V0 = V0
        inputs%pot_center = pot_center
        inputs%pot_width = pot_width
        
        inputs%max_iter = max_iter
        inputs%density_tol = density_tol
        inputs%energy_tol = energy_tol
        inputs%mixing_alpha = mixing_alpha
        inputs%verbose = verbose
        inputs%store_history = store_history
        
        inputs%output_prefix = output_prefix
        inputs%save_density = save_density
        inputs%save_eigenvalues = save_eigenvalues
        inputs%save_wavefunction = save_wavefunction
        
    end subroutine read_namelist_file

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
        
        if (inputs%Nup + inputs%Ndown > inputs%L) then
            print *, "ERROR: N = Nup + Ndown cannot exceed L"
            print *, "  N =", inputs%Nup + inputs%Ndown, ", L =", inputs%L
            ierr = ERROR_INVALID_INPUT
            return
        end if
        
        if (inputs%U < 0.0_dp) then
            print *, "ERROR: U must be non-negative, got U =", inputs%U
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
        
        ! Validate SCF parameters
        if (inputs%max_iter <= 0) then
            print *, "ERROR: max_iter must be positive"
            ierr = ERROR_INVALID_INPUT
            return
        end if
        
        if (inputs%mixing_alpha <= 0.0_dp .or. inputs%mixing_alpha > 1.0_dp) then
            print *, "ERROR: mixing_alpha must be in (0, 1]"
            ierr = ERROR_INVALID_INPUT
            return
        end if
        
    end subroutine validate_inputs

    !> Convert input_params_t to system_params_t
    !!
    !! @param[in]  inputs     Input parameters
    !! @param[out] sys_params System parameters
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
        sys_params%phase = inputs%phase
        
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
        scf_params%mixing_alpha = inputs%mixing_alpha
        scf_params%verbose = inputs%verbose
        scf_params%store_history = inputs%store_history
        
    end subroutine convert_to_scf_params

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
        print '(A)', "  --phase <real>    Twist angle for twisted BC (default: 0.0)"
        print '(A)', ""
        print '(A)', "Potential:"
        print '(A)', "  --potential <type>  Type: uniform/harmonic/barrier/... (default: uniform)"
        print '(A)', "  --V0 <real>         Potential strength (default: 0.0)"
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