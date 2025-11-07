module continuation
    use lsda_constants, only: dp, TWOPI, U_SMALL
    use bethe_equations, only: compute_energy
    use nonlinear_solvers, only: solve_newton
    implicit none
    private

    public :: estimate_dxdU
    public :: sweep_U_forward
    public :: sweep_U_backward
    public :: sweep_U_bidirectional

contains

    function estimate_dxdU(x_old, x_current, dU) result(dxdU)
        real(dp), intent(in) :: x_old(:), x_current(:), dU
        real(dp) :: dxdU(size(x_current))

        dxdU = (x_current - x_old) / dU
    end function estimate_dxdU

    subroutine sweep_U_forward(I, J, L, U_values, solutions, energies, converged_flags)
        real(dp), intent(in) :: U_values(:), I(:), J(:)
        integer, intent(in) :: L
        real(dp), intent(out) :: solutions(:, :), energies(:)
        logical, intent(out) :: converged_flags(:)

        real(dp), allocatable :: x(:), x_guess(:), dxdU(:), k(:)
        real(dp) :: U, dU, dU_old
        integer :: Nup, M, site, n_points
        logical :: converged

        Nup = size(I)
        M = size(J)
        n_points = size(U_values)

        allocate(x(Nup + M))
        allocate(x_guess(Nup + M))
        allocate(dxdU(Nup + M))
        allocate(k(Nup))

        ! For first site: site = 1
        U = U_values(1)

        if (U < U_SMALL) then
            ! Analytical solution for U = 0
            x_guess(1:Nup) = (TWOPI * I) / real(L, dp)
            x_guess(Nup + 1:Nup + M) = 0.0_dp
            x = x_guess
            converged = .true.
        else
            ! Initial guess: Fermi gas
            x_guess(1:Nup) = (TWOPI * I) / real(L, dp)
            x_guess(Nup + 1:Nup + M) = 0.0_dp

            x = x_guess
            call solve_newton(x, I, J, L, U, converged)
        end if

        solutions(:, 1) = x
        k = x(1:Nup)
        energies(1) = compute_energy(k)
        converged_flags(1) = converged

        ! For second site: site = 2
        U = U_values(2)
        x_guess = solutions(:, 1)
        x = x_guess
        call solve_newton(x, I, J, L, U, converged)

        solutions(:, 2) = x
        k = x(1:Nup)
        energies(2) = compute_energy(k)
        converged_flags(2) = converged

        dU_old = U_values(2) - U_values(1)
        do site = 3, n_points
            U = U_values(site)
            dU = U_values(site) - U_values(site - 1)

            ! PREDITOR
            dxdU = estimate_dxdU(solutions(:, site - 2), solutions(:, site - 1), dU_old)
            x_guess = solutions(:, site - 1) + dU * dxdU

            ! CORRETOR
            x = x_guess
            call solve_newton(x, I, J, L, U, converged)

            solutions(:, site) = x
            k = x(1:Nup)
            energies(site) = compute_energy(k)
            converged_flags(site) = converged

            if (.not. converged) then
                print *, "Warning: Failed at U =", U
            end if

            dU_old = dU
        end do

        deallocate(x, x_guess, dxdU, k)
    end subroutine

    subroutine sweep_U_backward(I, J, L, U_values, solutions, energies, converged_flags)
        real(dp), intent(in) :: U_values(:), I(:), J(:)
        integer, intent(in) :: L
        real(dp), intent(out) :: solutions(:, :), energies(:)
        logical, intent(out) :: converged_flags(:)

        real(dp), allocatable :: x(:), x_guess(:), dxdU(:), k(:)
        real(dp) :: U, dU, dU_old
        integer :: Nup, M, site, n_points
        logical :: converged

        Nup = size(I)
        M = size(J)
        n_points = size(U_values)

        allocate(x(Nup + M))
        allocate(x_guess(Nup + M))
        allocate(dxdU(Nup + M))
        allocate(k(Nup))

        ! For first site: site = n_points
        U = U_values(n_points)

        if (U < U_SMALL) then
            ! Analytical solution for U = 0
            x_guess(1:Nup) = (TWOPI * I) / real(L, dp)
            x_guess(Nup + 1:Nup + M) = 0.0_dp
            x = x_guess
            converged = .true.
        else
            ! Initial guess: Fermi gas
            x_guess(1:Nup) = (TWOPI * I) / real(L, dp)
            x_guess(Nup + 1:Nup + M) = 0.0_dp

            x = x_guess
            call solve_newton(x, I, J, L, U, converged)
        end if

        solutions(:, n_points) = x
        k = x(1:Nup)
        energies(n_points) = compute_energy(k)
        converged_flags(n_points) = converged

        ! For second site: site = n_points - 1
        U = U_values(n_points - 1)
        x_guess = solutions(:, n_points)
        x = x_guess
        call solve_newton(x, I, J, L, U, converged)

        solutions(:, n_points - 1) = x
        k = x(1:Nup)
        energies(n_points - 1) = compute_energy(k)
        converged_flags(n_points - 1) = converged

        dU_old = U_values(n_points) - U_values(n_points - 1)
        do site = n_points - 2, 1, -1
            U = U_values(site)
            dU = U_values(site + 1) - U_values(site)

            ! PREDITOR
            dxdU = estimate_dxdU(solutions(:, site + 2), solutions(:, site + 1), dU_old)
            x_guess = solutions(:, site + 1) - dU * dxdU

            ! CORRETOR
            x = x_guess
            call solve_newton(x, I, J, L, U, converged)

            solutions(:, site) = x
            k = x(1:Nup)
            energies(site) = compute_energy(k)
            converged_flags(site) = converged

            if (.not. converged) then
                print *, "Warning: Failed at U =", U
            end if

            dU_old = dU
        end do

        deallocate(x, x_guess, dxdU, k)
    end subroutine

    subroutine sweep_U_bidirectional(I, J, L, U_values, solutions, energies, converged_flags)
        real(dp), intent(in) :: U_values(:), I(:), J(:)
        integer, intent(in) :: L
        real(dp), intent(out) :: solutions(:, :), energies(:)
        logical, intent(out) :: converged_flags(:)

        real(dp), allocatable :: sol_fwd(:, :), sol_bwd(:, :), E_fwd(:), E_bwd(:)
        logical, allocatable :: flags_fwd(:), flags_bwd(:)
        real(dp) :: max_diff
        integer :: Nup, M, n_points

        Nup = size(I)
        M = size(J)
        n_points = size(U_values)

        allocate(sol_fwd(Nup+M, n_points))
        allocate(E_fwd(n_points))
        allocate(flags_fwd(n_points))
        allocate(sol_bwd(Nup+M, n_points))
        allocate(E_bwd(n_points))
        allocate(flags_bwd(n_points))

        ! 1. Forward sweep
        call sweep_U_forward(I, J, L, U_values, sol_fwd, E_fwd, flags_fwd)
        
        ! 2. Backward sweep
        call sweep_U_backward(I, J, L, U_values, sol_bwd, E_bwd, flags_bwd)

        ! 3. Average (refinment)
        solutions = (sol_fwd + sol_bwd) / 2.0_dp
        energies = (E_fwd + E_bwd) / 2.0_dp

        ! 4. Convergence: both must have converged
        converged_flags = flags_fwd .and. flags_bwd
        
        ! 5. Check consistency
        max_diff = maxval(abs(E_fwd - E_bwd))
        if (max_diff > 1.0e-6_dp) then
            print *, "Warning: Forward/backward differ by", max_diff
        end if
    end subroutine
end module continuation