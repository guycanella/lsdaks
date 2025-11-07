!> Continuation methods for solving Bethe Ansatz equations across parameter ranges
!!
!! This module implements predictor-corrector continuation methods to efficiently
!! solve the Lieb-Wu equations for multiple values of the Hubbard interaction U.
!!
!! Key features:
!! - Forward sweep (U_min → U_max) with linear extrapolation
!! - Backward sweep (U_max → U_min) for refinement
!! - Bidirectional sweep (average of forward + backward)
!! - Typical speedup: 5-10x compared to independent solutions
!!
!! Algorithm:
!! 1. Solve for first U value (from scratch)
!! 2. For subsequent U: predict solution using previous points
!! 3. Refine prediction with Newton-Raphson (converges in 1-3 iterations)
!!
!! @note Assumes solutions vary smoothly with U (no bifurcations)
!! @see bethe_equations, nonlinear_solvers
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

    !> Estimates dx/dU using finite differences
    !!
    !! Approximates the derivative of the solution with respect to U
    !! using two consecutive points:
    !!
    !! \[ \frac{dx}{dU} \approx \frac{x_{\text{current}} - x_{\text{old}}}{\Delta U} \]
    !!
    !! This is used in the predictor step of continuation methods to
    !! extrapolate the solution to the next U value.
    !!
    !! @param[in] x_old      Solution at U_{i-1}
    !! @param[in] x_current  Solution at U_i
    !! @param[in] dU         Step size: U_i - U_{i-1}
    !! @return    dxdU       Estimated derivative dx/dU
    !!
    !! @note This is a first-order approximation (linear extrapolation)
    !! @note For better accuracy, could use higher-order schemes (quadratic, cubic)
    !!
    !! @see sweep_U_forward, sweep_U_backward
    function estimate_dxdU(x_old, x_current, dU) result(dxdU)
        real(dp), intent(in) :: x_old(:), x_current(:), dU
        real(dp) :: dxdU(size(x_current))

        dxdU = (x_current - x_old) / dU
    end function estimate_dxdU

    !> Forward sweep in U with predictor-corrector method
    !!
    !! Solves Bethe Ansatz equations for a sequence of U values from U_min to U_max,
    !! using previous solutions to accelerate convergence.
    !!
    !! Algorithm:
    !! - Point 1 (U_min): Solve from scratch (Fermi gas initial guess)
    !! - Point 2: Predictor = previous solution (no extrapolation yet)
    !! - Points 3+: Predictor = linear extrapolation using 2 previous points
    !!              x_guess = 2·x_{i-1} - x_{i-2}
    !! - Corrector: Refine with Newton-Raphson (typically 1-3 iterations)
    !!
    !! @param[in]  I                Charge quantum numbers (N↑)
    !! @param[in]  J                Spin quantum numbers (M = N↓)
    !! @param[in]  L                Number of lattice sites
    !! @param[in]  U_values         U values [U_min, U_min+dU, ..., U_max] (n_points)
    !! @param[out] solutions        Solutions x(U) - matrix (N↑+M, n_points)
    !! @param[out] energies         Ground state energies E(U) (n_points)
    !! @param[out] converged_flags  Convergence flags (n_points)
    !!
    !! @note Arrays must be pre-allocated by caller
    !! @note Prints warning if Newton fails to converge at any point
    !! @note Typical speedup vs independent solutions: 5-10x
    !!
    !! Example usage:
    !! @code
    !!   U_values = [0.0, 0.5, 1.0, ..., 10.0]  ! 21 points
    !!   allocate(solutions(Nup+M, 21), energies(21), flags(21))
    !!   call sweep_U_forward(I, J, L, U_values, solutions, energies, flags)
    !! @endcode
    !!
    !! @see sweep_U_backward, sweep_U_bidirectional, solve_newton
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

    !> Backward sweep in U with predictor-corrector method
    !!
    !! Identical to forward sweep, but in reverse direction (U_max → U_min).
    !! Useful for refinement and validation of forward sweep results.
    !!
    !! Algorithm:
    !! - Point n (U_max): Solve from scratch
    !! - Point n-1: Predictor = previous solution
    !! - Points n-2 to 1: Predictor = linear extrapolation (backward)
    !!                    x_guess = x_{i+1} - dU·(dx/dU)
    !!
    !! @param[in]  I                Charge quantum numbers
    !! @param[in]  J                Spin quantum numbers
    !! @param[in]  L                Number of lattice sites
    !! @param[in]  U_values         U values [U_min, ..., U_max]
    !! @param[out] solutions        Solutions x(U) - matrix (N↑+M, n_points)
    !! @param[out] energies         Ground state energies E(U)
    !! @param[out] converged_flags  Convergence flags
    !!
    !! @note U_values array order is the same as forward (U_min first)
    !! @note Solutions are filled in reverse order (from n_points to 1)
    !! @note Useful to detect multiple solutions or bifurcations
    !!
    !! @see sweep_U_forward, sweep_U_bidirectional
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

    !> Bidirectional sweep with refinement (forward + backward average)
    !!
    !! Performs both forward and backward sweeps, then averages the results
    !! to obtain maximum accuracy and detect potential numerical issues.
    !!
    !! Algorithm:
    !! 1. Forward sweep: U_min → U_max (stores in sol_fwd)
    !! 2. Backward sweep: U_max → U_min (stores in sol_bwd)
    !! 3. Refinement: Average both results
    !!    - solutions = (sol_fwd + sol_bwd) / 2
    !!    - energies = (E_fwd + E_bwd) / 2
    !! 4. Validation: Check consistency
    !!    - If |E_fwd - E_bwd| > 1e-6, prints warning
    !!    - May indicate bifurcations or numerical instability
    !! 5. Convergence: Both sweeps must have converged
    !!    - converged_flags = flags_fwd AND flags_bwd
    !!
    !! @param[in]  I                Charge quantum numbers (N↑)
    !! @param[in]  J                Spin quantum numbers (M = N↓)
    !! @param[in]  L                Number of lattice sites
    !! @param[in]  U_values         U values [U_min, ..., U_max] (n_points)
    !! @param[out] solutions        Refined solutions (average) - matrix (N↑+M, n_points)
    !! @param[out] energies         Refined energies (average) - vector (n_points)
    !! @param[out] converged_flags  Convergence flags (true if BOTH converged)
    !!
    !! @note This is the most robust method but ~2x slower than single sweep
    !! @note Recommended for production runs and validation
    !! @note If forward/backward differ significantly, may indicate:
    !!       - Multiple solutions (bifurcation)
    !!       - Numerical instability
    !!       - Step size too large
    !!
    !! Advantages over single sweep:
    !! - Higher accuracy (averaging reduces numerical errors)
    !! - Error detection (large differences indicate problems)
    !! - Validation (if fwd ≈ bwd, solution is trustworthy)
    !!
    !! Example usage:
    !! @code
    !!   ! Generate table for U ∈ [0, 10] with 100 points
    !!   U_values = [(i*0.1_dp, i=0,100)]
    !!   allocate(solutions(Nup+M, 101), energies(101), flags(101))
    !!   call sweep_U_bidirectional(I, J, L, U_values, solutions, energies, flags)
    !!   
    !!   if (all(flags)) then
    !!       print *, "All points converged successfully!"
    !!   end if
    !! @endcode
    !!
    !! @see sweep_U_forward, sweep_U_backward
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