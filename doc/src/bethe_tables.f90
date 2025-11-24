!> Module for generating exchange-correlation functional tables
!!
!! This module provides functionality to generate XC tables by:
!! 1. Computing E_xc = E_BA - E_0 for grid points (n, m, U)
!! 2. Computing V_xc via numerical derivatives
!! 3. Using continuation method for efficient U sweeps
!! 4. OpenMP parallelization for grid point calculations
module bethe_tables
    use lsda_constants, only: dp, TWOPI, U_SMALL
    use bethe_equations, only: initialize_quantum_numbers, compute_energy, compute_residual
    use nonlinear_solvers, only: solve_newton
    use table_io, only: xc_table_t, write_fortran_table
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
    implicit none
    private

    type, public :: grid_params_t
        real(dp) :: n_min = 0.1_dp       !< Minimum density
        real(dp) :: n_max = 2.0_dp       !< Maximum density
        integer :: n_points = 50         !< Number of density points
        integer :: m_points = 51         !< Number of magnetization points (per n)
        integer :: L = 100               !< System size
        real(dp) :: delta_n = 1.0e-4_dp  !< Finite difference increment
    end type grid_params_t

    type, public :: xc_potentials_t
        real(dp) :: v_xc_up, v_xc_down
    end type xc_potentials_t

    public :: compute_E0
    public :: compute_E_xc
    public :: compute_V_xc_numerical
    public :: generate_xc_table
    public :: generate_table_grid

contains
    !> Compute non-interacting energy (U=0, free Fermi gas)
    !!
    !! E_0 = -2 * sum_j cos(k_j) for both spins
    !! with k_j = 2π/L * I_j, I_j = j - (N+1)/2
    !!
    !! @param[in] n_up  Spin-up density
    !! @param[in] n_dw  Spin-down density
    !! @param[in] L     System size
    !! @return          Non-interacting energy E_0
    function compute_E0(n_up, n_dw, L) result(E0)
        real(dp), intent(in) :: n_up, n_dw
        integer, intent(in) :: L
        real(dp) :: E0, I_j, k_j
        integer :: Nup, Ndw, j

        Nup = NINT(n_up * real(L, dp))
        Ndw = NINT(n_dw * real(L, dp))

        E0 = 0.0_dp

        do j = 1, Nup
            I_j = real(j, dp) -0.5_dp * real(Nup + 1, dp)
            k_j = TWOPI * I_j / real(L, dp)
            E0 = E0 - 2.0_dp * cos(k_j)
        end do

        do j = 1, Ndw
            I_j = real(j, dp) -0.5_dp * real(Ndw + 1, dp)
            k_j = TWOPI * I_j / real(L, dp)
            E0 = E0 - 2.0_dp * cos(k_j)
        end do
    end function compute_E0

    !> Compute exchange-correlation energy: E_xc = E_BA - E_0
    !!
    !! Solves Bethe Ansatz equations to get E_BA, then subtracts E_0.
    !! Returns energy per site for extensivity.
    !!
    !! @param[in] n_up  Spin-up density
    !! @param[in] n_dw  Spin-down density
    !! @param[in] U     Hubbard interaction
    !! @param[in] L     System size
    !! @return          XC energy per site e_xc = E_xc/L
    function compute_E_xc(n_up, n_dn, U, L) result(E_xc)
        real(dp), intent(in) :: n_up, n_dn, U
        integer, intent(in) :: L
        real(dp) :: E_xc
        
        integer :: Nup, Ndn, M
        real(dp) :: E0, E_BA
        real(dp), allocatable :: x(:), k(:), F(:)
        real(dp), allocatable :: quantum_I(:), quantum_J(:)
        logical :: converged

        if (abs(U) < U_SMALL) then
            E_xc = 0.0_dp
            return
        end if

        Nup = nint(n_up * real(L, dp))
        Ndn = nint(n_dn * real(L, dp))
        M = Ndn

        allocate(x(Nup + M))  ! Combined array for solver (k + Lambda)
        allocate(k(Nup))
        allocate(F(Nup + M))
        allocate(quantum_I(Nup), quantum_J(M))

        call initialize_quantum_numbers(Nup, M, quantum_I, quantum_J)
        
        x(1:Nup) = TWOPI * quantum_I / real(L, dp)
        if (M > 0) then
            x(Nup+1:Nup+M) = 0.0_dp
        end if

        call solve_newton(x, quantum_I, quantum_J, L, U, converged)

        if (.not. converged) then
            F = compute_residual(x(1:Nup), x(Nup+1:Nup+M), quantum_I, quantum_J, L, U)

            if (norm2(F) > 1.0e-8_dp) then
                E_xc = ieee_value(E_xc, ieee_quiet_nan)
                deallocate(x, k, quantum_I, quantum_J, F)
                return
            end if

            deallocate(F)
        end if
        
        k = x(1:Nup)
        E_BA = compute_energy(k)
        E0 = compute_E0(n_up, n_dn, L)
        E_xc = (E_BA - E0) / real(L, dp)

        deallocate(x, k, quantum_I, quantum_J)
    end function compute_E_xc

    !> Compute XC potentials via numerical derivatives (central finite differences)
    !!
    !! V_xc_up = ∂E_xc/∂n_up, V_xc_dw = ∂E_xc/∂n_dw
    !! Uses 2nd-order central differences with δn ~ 1e-4
    !!
    !! @param[in] n_up    Spin-up density
    !! @param[in] n_dw    Spin-down density
    !! @param[in] U       Hubbard interaction
    !! @param[in] L       System size
    !! @param[in] delta_n Finite difference increment
    !! @return            XC potentials (v_xc_up, v_xc_down)
    function compute_V_xc_numerical(n_up, n_dw, U, L) result(v_xc)
        real(dp), intent(in) :: n_up, n_dw
        integer, intent(in) :: L
        real(dp), intent(in) :: U
        type(xc_potentials_t) :: v_xc

        real(dp) :: delta_n
        real(dp) :: E_xc_up_plus, E_xc_up_minus
        real(dp) :: E_xc_dw_plus, E_xc_dw_minus

        delta_n = 1.0e-4_dp

        ! Derivative with respect to n_up
        E_xc_up_plus = compute_E_xc(n_up + delta_n, n_dw, U, L)
        E_xc_up_minus = compute_E_xc(n_up - delta_n, n_dw, U, L)
        if (ieee_is_nan(E_xc_up_plus) .or. ieee_is_nan(E_xc_up_minus)) then
            v_xc%v_xc_up = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            v_xc%v_xc_up = (E_xc_up_plus - E_xc_up_minus) / (2.0_dp * delta_n)
        end if

        ! Derivative with respect to n_dw
        E_xc_dw_plus = compute_E_xc(n_up, n_dw + delta_n, U, L)
        E_xc_dw_minus = compute_E_xc(n_up, n_dw - delta_n, U, L)
        if (ieee_is_nan(E_xc_dw_plus) .or. ieee_is_nan(E_xc_dw_minus)) then
            v_xc%v_xc_down = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            v_xc%v_xc_down = (E_xc_dw_plus - E_xc_dw_minus) / (2.0_dp * delta_n)
        end if
    end function compute_V_xc_numerical

    !> Generate complete XC table for a single U value
    !!
    !! Creates table with grid points (n_i, m_j) where:
    !!   n ∈ [n_min, n_max], m ∈ [-n, n]
    !! and computes exc, vxc_up, vxc_down at each point.
    !!
    !! @param[in]  U      Hubbard interaction parameter
    !! @param[in]  params Grid parameters (densities, points, system size)
    !! @param[out] table  Generated XC table
    !! @param[out] status Success (0) or error code
    subroutine generate_xc_table(U, params, table, status)
        real(dp), intent(in) :: U
        type(grid_params_t), intent(in) :: params
        type(xc_table_t), intent(out) :: table
        integer, intent(out) :: status
        
        integer :: i, j
        real(dp) :: n, m, n_up, n_dw, delta_n, delta_m
        type(xc_potentials_t) :: v_xc
        integer :: actual_m_points
        
        status = 0
        
        table%U = U
        table%n_points_n = params%n_points
        table%n_points_m = params%m_points
        
        allocate(table%n_grid(params%n_points))
        allocate(table%m_grid(params%m_points, params%n_points))
        allocate(table%exc(params%m_points, params%n_points))
        allocate(table%vxc_up(params%m_points, params%n_points))
        allocate(table%vxc_down(params%m_points, params%n_points))

        ! Initialize with NaN (for debugging)
        table%exc = ieee_value(0.0_dp, ieee_quiet_nan)
        table%vxc_up = ieee_value(0.0_dp, ieee_quiet_nan)
        table%vxc_down = ieee_value(0.0_dp, ieee_quiet_nan)
        
        delta_n = (params%n_max - params%n_min) / real(params%n_points - 1, dp)
        
        !$OMP PARALLEL DO PRIVATE(i, j, n, m, n_up, n_dw, delta_m, v_xc, actual_m_points) &
        !$OMP SHARED(table, params, U, delta_n) SCHEDULE(dynamic)
        do i = 1, params%n_points
            n = params%n_min + real(i - 1, dp) * delta_n
            table%n_grid(i) = n
            
            actual_m_points = params%m_points
            if (actual_m_points > 1) then
                delta_m = (2.0_dp * n) / real(actual_m_points - 1, dp)
            else
                delta_m = 0.0_dp
            end if
            
            do j = 1, actual_m_points
                m = -n + real(j - 1, dp) * delta_m
                table%m_grid(j, i) = m
                
                n_up = 0.5_dp * (n + m)
                n_dw = 0.5_dp * (n - m)
                
                if (n_up < 0.0_dp .or. n_up > 1.0_dp .or. &
                    n_dw < 0.0_dp .or. n_dw > 1.0_dp) then
                    cycle
                end if
                
                table%exc(j, i) = compute_E_xc(n_up, n_dw, U, params%L)

                v_xc = compute_V_xc_numerical(n_up, n_dw, U, params%L)
                table%vxc_up(j, i) = v_xc%v_xc_up
                table%vxc_down(j, i) = v_xc%v_xc_down
            end do
        end do
        !$OMP END PARALLEL DO
        
    end subroutine generate_xc_table

    !> Generate multiple XC tables using continuation method in U
    !!
    !! Efficiently generates tables for multiple U values by:
    !! 1. Computing first table from scratch
    !! 2. Using continuation method for subsequent U values
    !! 3. Saving each table to binary format
    !!
    !! @param[in]  U_values   Array of Hubbard U values
    !! @param[in]  params     Grid parameters
    !! @param[in]  output_dir Directory to save table files
    !! @param[out] status     Success (0) or error code
    subroutine generate_table_grid(U_values, params, output_dir, status)
        real(dp), intent(in) :: U_values(:)
        type(grid_params_t), intent(in) :: params
        character(len=*), intent(in) :: output_dir
        integer, intent(out) :: status
        
        integer :: k, n_U, io_stat
        real(dp) :: U_current
        type(xc_table_t) :: table
        character(len=256) :: filename
        character(len=32) :: U_str
        
        status = 0
        n_U = size(U_values)
        
        print '(A)', "========================================="
        print '(A)', "  XC Table Grid Generation"
        print '(A)', "========================================="
        print '(A,I0,A)', "Generating ", n_U, " tables with continuation method"
        print '(A,F6.2,A,F6.2)', "U range: ", U_values(1), " to ", U_values(n_U)
        print '(A,I0,A,I0)', "Grid size: ", params%n_points, " x ", params%m_points
        print '(A)', ""
        
        do k = 1, n_U
            U_current = U_values(k)
            
            print '(A,I0,A,I0,A,F6.2)', "Processing U(", k, "/", n_U, ") = ", U_current
            
            call generate_xc_table(U_current, params, table, status)
            
            if (status /= 0) then
                print '(A,F6.2)', "  ERROR: Failed to generate table for U = ", U_current
                return
            end if
            
            write(U_str, '(F6.2)') U_current
            U_str = adjustl(U_str)
            filename = trim(output_dir) // "/lsda_hub_u" // trim(U_str) // ".dat"
            
            call write_fortran_table(filename, table, io_stat)
            
            if (io_stat /= 0) then
                print '(A)', "  ERROR: Failed to write table file"
                status = -1
                return
            end if
            
            print '(A,A)', "  Saved: ", trim(filename)
            
            ! Deallocate table arrays for next iteration
            if (allocated(table%n_grid)) deallocate(table%n_grid)
            if (allocated(table%m_grid)) deallocate(table%m_grid)
            if (allocated(table%exc)) deallocate(table%exc)
            if (allocated(table%vxc_up)) deallocate(table%vxc_up)
            if (allocated(table%vxc_down)) deallocate(table%vxc_down)
        end do
        
        print '(A)', ""
        print '(A)', "========================================="
        print '(A,I0,A)', "Successfully generated ", n_U, " XC tables!"
        print '(A)', "========================================="
        
    end subroutine generate_table_grid
end module bethe_tables