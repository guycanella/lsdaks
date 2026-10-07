!> Continuation methods for solving Bethe Ansatz equations across parameter ranges
!!
!! @warning **Finite-size validation tool, not production.** Since phase 4.5
!!          the XC table generator uses the thermodynamic-limit integral
!!          equations of `lieb_wu_integral`, so this module and the discrete
!!          solvers it drives (`bethe_equations`, `nonlinear_solvers`) are no
!!          longer on the production path.  They are kept to cross-check the
!!          integral-equation results against explicit finite-L solutions.
!!          The parity of the quantum numbers and the degenerate initial guess
!!          for the spin rapidities are still known to be wrong for part of
!!          the parameter range, which is why every `solve_newton` return is
!!          validated here before the solution is used.
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
    use lsda_errors, only: ERROR_SUCCESS
    use bethe_equations, only: compute_energy
    use nonlinear_solvers, only: solve_newton
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
    implicit none
    private

    public :: estimate_dxdU
    public :: sweep_U_forward
    public :: sweep_U_backward
    public :: sweep_U_bidirectional
    public :: merge_sweeps

    !> Largest forward/backward energy difference accepted without a warning.
    real(dp), parameter, public :: E_CROSS_TOL = 1.0e-6_dp

contains

    !> Store one sweep point, rejecting an unconverged Newton return.
    !!
    !! `solve_newton` can return a vector that is merely the last iterate: it
    !! sets `converged = .false.`, or a non-zero `solver_status` when the
    !! linear algebra itself failed.  Storing such a vector poisons the energy,
    !! the reported flag, and - worst of all - the predictor of every later
    !! point.  This routine therefore writes NaN instead, and reports whether
    !! the point may be used as a predictor.
    !!
    !! @param[in]  x              Newton return vector
    !! @param[in]  Nup            Number of charge rapidities in `x`
    !! @param[in]  converged      Convergence flag from `solve_newton`
    !! @param[in]  solver_status  Error code from `solve_newton`
    !! @param[out] solution       Stored solution (NaN when rejected)
    !! @param[out] energy         Stored energy (NaN when rejected)
    !! @param[out] accepted       `.true.` when the point may be reused
    subroutine store_point(x, Nup, converged, solver_status, solution, energy, accepted)
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: Nup, solver_status
        logical, intent(in) :: converged
        real(dp), intent(out) :: solution(:), energy
        logical, intent(out) :: accepted

        accepted = converged .and. (solver_status == ERROR_SUCCESS)

        if (accepted) then
            solution = x
            energy = compute_energy(x(1:Nup))
        else
            solution = ieee_value(0.0_dp, ieee_quiet_nan)
            energy = ieee_value(0.0_dp, ieee_quiet_nan)
        end if
    end subroutine store_point

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

        real(dp), allocatable :: x(:), x_guess(:), dxdU(:), x_fermi(:), x_good(:)
        real(dp) :: U, dU, dU_old
        integer :: Nup, M, site, n_points, solver_status
        logical :: converged, accepted, have_good, have_two

        Nup = size(I)
        M = size(J)
        n_points = size(U_values)

        allocate(x(Nup + M))
        allocate(x_guess(Nup + M))
        allocate(dxdU(Nup + M))
        allocate(x_fermi(Nup + M))
        allocate(x_good(Nup + M))

        x_fermi(1:Nup) = (TWOPI * I) / real(L, dp)
        x_fermi(Nup + 1:Nup + M) = 0.0_dp
        x_good = x_fermi
        have_good = .false.
        have_two = .false.
        solver_status = ERROR_SUCCESS

        ! For first site: site = 1
        U = U_values(1)

        x = x_fermi
        if (U < U_SMALL) then
            ! Analytical solution for U = 0
            converged = .true.
        else
            call solve_newton(x, I, J, L, U, converged, solver_status)
        end if

        call store_point(x, Nup, converged, solver_status, solutions(:, 1), &
                         energies(1), accepted)
        converged_flags(1) = accepted
        if (accepted) then
            x_good = x
            have_good = .true.
        else
            print *, "Warning: Failed at U =", U
        end if

        ! For second site: site = 2
        U = U_values(2)
        x = x_good
        call solve_newton(x, I, J, L, U, converged, solver_status)

        call store_point(x, Nup, converged, solver_status, solutions(:, 2), &
                         energies(2), accepted)
        converged_flags(2) = accepted
        if (accepted) then
            have_two = have_good
            x_good = x
            have_good = .true.
        else
            have_two = .false.
            print *, "Warning: Failed at U =", U
        end if

        dU_old = U_values(2) - U_values(1)
        do site = 3, n_points
            U = U_values(site)
            dU = U_values(site) - U_values(site - 1)

            ! PREDITOR - only extrapolate through points that were accepted;
            ! a rejected point holds NaN and would destroy the guess.
            if (have_two .and. converged_flags(site - 1) .and. converged_flags(site - 2)) then
                dxdU = estimate_dxdU(solutions(:, site - 2), solutions(:, site - 1), dU_old)
                x_guess = solutions(:, site - 1) + dU * dxdU
            else if (have_good) then
                x_guess = x_good
            else
                x_guess = x_fermi
            end if

            ! CORRETOR
            x = x_guess
            call solve_newton(x, I, J, L, U, converged, solver_status)

            call store_point(x, Nup, converged, solver_status, solutions(:, site), &
                             energies(site), accepted)
            converged_flags(site) = accepted
            if (accepted) then
                have_two = have_good
                x_good = x
                have_good = .true.
            else
                have_two = .false.
                print *, "Warning: Failed at U =", U
            end if

            dU_old = dU
        end do

        deallocate(x, x_guess, dxdU, x_fermi, x_good)
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

        real(dp), allocatable :: x(:), x_guess(:), dxdU(:), x_fermi(:), x_good(:)
        real(dp) :: U, dU, dU_old
        integer :: Nup, M, site, n_points, solver_status
        logical :: converged, accepted, have_good, have_two

        Nup = size(I)
        M = size(J)
        n_points = size(U_values)

        allocate(x(Nup + M))
        allocate(x_guess(Nup + M))
        allocate(dxdU(Nup + M))
        allocate(x_fermi(Nup + M))
        allocate(x_good(Nup + M))

        x_fermi(1:Nup) = (TWOPI * I) / real(L, dp)
        x_fermi(Nup + 1:Nup + M) = 0.0_dp
        x_good = x_fermi
        have_good = .false.
        have_two = .false.
        solver_status = ERROR_SUCCESS

        ! For first site: site = n_points
        U = U_values(n_points)

        x = x_fermi
        if (U < U_SMALL) then
            ! Analytical solution for U = 0
            converged = .true.
        else
            call solve_newton(x, I, J, L, U, converged, solver_status)
        end if

        call store_point(x, Nup, converged, solver_status, solutions(:, n_points), &
                         energies(n_points), accepted)
        converged_flags(n_points) = accepted
        if (accepted) then
            x_good = x
            have_good = .true.
        else
            print *, "Warning: Failed at U =", U
        end if

        ! For second site: site = n_points - 1
        U = U_values(n_points - 1)
        x = x_good
        call solve_newton(x, I, J, L, U, converged, solver_status)

        call store_point(x, Nup, converged, solver_status, solutions(:, n_points - 1), &
                         energies(n_points - 1), accepted)
        converged_flags(n_points - 1) = accepted
        if (accepted) then
            have_two = have_good
            x_good = x
            have_good = .true.
        else
            have_two = .false.
            print *, "Warning: Failed at U =", U
        end if

        dU_old = U_values(n_points) - U_values(n_points - 1)
        do site = n_points - 2, 1, -1
            U = U_values(site)
            dU = U_values(site + 1) - U_values(site)

            ! PREDITOR - see the forward sweep: a rejected point holds NaN.
            if (have_two .and. converged_flags(site + 1) .and. converged_flags(site + 2)) then
                ! The arguments are ordered by increasing U so that dxdU is
                ! the physical dx/dU; the subtraction below then extrapolates
                ! from U(site+1) towards the lower-U point.
                dxdU = estimate_dxdU(solutions(:, site + 1), solutions(:, site + 2), dU_old)
                x_guess = solutions(:, site + 1) - dU * dxdU
            else if (have_good) then
                x_guess = x_good
            else
                x_guess = x_fermi
            end if

            ! CORRETOR
            x = x_guess
            call solve_newton(x, I, J, L, U, converged, solver_status)

            call store_point(x, Nup, converged, solver_status, solutions(:, site), &
                             energies(site), accepted)
            converged_flags(site) = accepted
            if (accepted) then
                have_two = have_good
                x_good = x
                have_good = .true.
            else
                have_two = .false.
                print *, "Warning: Failed at U =", U
            end if

            dU_old = dU
        end do

        deallocate(x, x_guess, dxdU, x_fermi, x_good)
    end subroutine

    !> Bidirectional sweep with refinement (forward + backward average)
    !!
    !! Performs both forward and backward sweeps, then averages the results
    !! to obtain maximum accuracy and detect potential numerical issues.
    !!
    !! Algorithm:
    !! 1. Forward sweep: U_min → U_max (stores in sol_fwd)
    !! 2. Backward sweep: U_max → U_min (stores in sol_bwd)
    !! 3. Refinement: Average the two results **where both were accepted**.
    !!    A rejected point holds NaN (see `store_point`), so an unconditional
    !!    average would turn a perfectly good forward point into NaN just
    !!    because the backward sweep failed at the same `U`.  Where only one
    !!    direction was accepted, that direction is used as is.
    !! 4. Validation: `max |E_fwd - E_bwd|` over the points accepted in both
    !!    directions; a difference above 1e-6 prints a warning and may
    !!    indicate bifurcations or numerical instability.  The comparison is
    !!    written so that a NaN which escapes the mask still fires the
    !!    warning, because `NaN > tol` is false in IEEE and would otherwise
    !!    silence the only cross-check this module has.
    !! 5. Convergence: `converged_flags = flags_fwd AND flags_bwd`; a point
    !!    accepted by a single direction is still reported as not converged,
    !!    because it has no cross-check behind it.
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
        logical :: inconsistent

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

        ! 3.-5. Merge, cross-check, report.
        call merge_sweeps(sol_fwd, E_fwd, flags_fwd, sol_bwd, E_bwd, flags_bwd, &
                          solutions, energies, converged_flags, max_diff, inconsistent)

        if (inconsistent) print *, "Warning: Forward/backward differ by", max_diff

        if (.not. all(converged_flags)) then
            print *, "Warning: bidirectional sweep has", &
                     count(.not. converged_flags), "point(s) without a cross-check"
        end if
    end subroutine

    !> Merge a forward and a backward sweep into one refined result.
    !!
    !! Split out of `sweep_U_bidirectional` so that the masking can be tested
    !! without having to drive `solve_newton` into failure.
    !!
    !! Per point:
    !! * accepted in both directions: arithmetic mean, and the point takes part
    !!   in the forward/backward cross-check;
    !! * accepted in one direction only: that direction is copied through, and
    !!   the point does **not** take part in the cross-check - it has no second
    !!   opinion behind it - but its good value is not destroyed either;
    !! * accepted in neither: NaN.
    !!
    !! `inconsistent` is computed with `.not. (x <= tol)` throughout, so a NaN
    !! that slipped past the mask reports as inconsistent instead of silently
    !! passing (`NaN > tol` is `.false.` in IEEE arithmetic).
    !!
    !! @param[in]  sol_fwd          Forward solutions
    !! @param[in]  E_fwd            Forward energies
    !! @param[in]  flags_fwd        Forward acceptance flags
    !! @param[in]  sol_bwd          Backward solutions
    !! @param[in]  E_bwd            Backward energies
    !! @param[in]  flags_bwd        Backward acceptance flags
    !! @param[out] solutions        Merged solutions
    !! @param[out] energies         Merged energies
    !! @param[out] converged_flags  `.true.` only where both directions agreed to run
    !! @param[out] max_diff         Largest `|E_fwd - E_bwd|` over cross-checked points
    !! @param[out] inconsistent     `.true.` when `max_diff` exceeds `E_CROSS_TOL`
    subroutine merge_sweeps(sol_fwd, E_fwd, flags_fwd, sol_bwd, E_bwd, flags_bwd, &
                            solutions, energies, converged_flags, max_diff, inconsistent)
        real(dp), intent(in) :: sol_fwd(:, :), E_fwd(:), sol_bwd(:, :), E_bwd(:)
        logical, intent(in) :: flags_fwd(:), flags_bwd(:)
        real(dp), intent(out) :: solutions(:, :), energies(:)
        logical, intent(out) :: converged_flags(:)
        real(dp), intent(out) :: max_diff
        logical, intent(out) :: inconsistent

        integer :: site
        real(dp) :: diff
        logical :: any_both, nan_seen

        max_diff = 0.0_dp
        any_both = .false.
        nan_seen = .false.

        do site = 1, size(E_fwd)
            if (flags_fwd(site) .and. flags_bwd(site)) then
                solutions(:, site) = 0.5_dp * (sol_fwd(:, site) + sol_bwd(:, site))
                energies(site) = 0.5_dp * (E_fwd(site) + E_bwd(site))
                diff = abs(E_fwd(site) - E_bwd(site))
                ! `max_diff` alone is not enough to carry a NaN: the next
                ! iteration would compare against it (`.not. (diff <= NaN)` is
                ! `.true.`) and overwrite it with a finite value, silencing the
                ! warning.  The NaN flag is therefore accumulated separately.
                nan_seen = nan_seen .or. ieee_is_nan(diff)
                if (.not. (diff <= max_diff)) max_diff = diff
                any_both = .true.
            else if (flags_fwd(site)) then
                solutions(:, site) = sol_fwd(:, site)
                energies(site) = E_fwd(site)
            else if (flags_bwd(site)) then
                solutions(:, site) = sol_bwd(:, site)
                energies(site) = E_bwd(site)
            else
                solutions(:, site) = ieee_value(0.0_dp, ieee_quiet_nan)
                energies(site) = ieee_value(0.0_dp, ieee_quiet_nan)
            end if
        end do

        converged_flags = flags_fwd .and. flags_bwd
        inconsistent = any_both .and. (nan_seen .or. .not. (max_diff <= E_CROSS_TOL))
    end subroutine merge_sweeps
end module continuation
