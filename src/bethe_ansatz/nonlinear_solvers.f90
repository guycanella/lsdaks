!> Nonlinear solvers for the Bethe Ansatz equations
!!
!! This module provides robust numerical methods to solve the nonlinear system
!! F(x) = 0 arising from the Lieb-Wu equations of the 1D Hubbard model.
!!
!! Main components:
!! - Newton-Raphson method with analytical Jacobian
!! - Armijo line search for global convergence
!! - LAPACK interface for linear systems
!!
!! @note All solvers use double precision (real64) arithmetic
module nonlinear_solvers
    use lsda_constants, only: dp, TWOPI, U_SMALL, NEWTON_TOL, NEWTON_MAX_ITER
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use bethe_equations
    implicit none
    private

    public :: solve_newton
    public :: solve_linear_system
    public :: line_search

    !> Armijo condition parameter (c ∈ (0,1))
    !! Smaller values = more aggressive line search
    real(dp), parameter :: ARMIJO_C = 1.0e-4_dp
    !> Maximum iterations for line search backtracking
    integer, parameter :: MAX_LS_ITER = 20

    !> LAPACK interface for solving linear systems A·x = b
    !! Uses LU decomposition with partial pivoting (DGESV)
    interface
        subroutine DGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
            import :: dp
            integer, intent(in) :: N, NRHS, LDA, LDB
            real(dp), intent(inout) :: A(LDA, *)
            integer, intent(out) :: IPIV(*)
            real(dp), intent(inout) :: B(LDB, *)
            integer, intent(out) :: INFO
        end subroutine DGESV
    end interface

contains

    !> Solves linear system A·x = b using LAPACK DGESV
    !!
    !! Wrapper for LAPACK's DGESV routine (LU decomposition with partial pivoting).
    !! Preserves input arrays by creating internal copies.
    !!
    !! @param[in]  A  Coefficient matrix (N×N, must be square and non-singular)
    !! @param[out] x  Solution vector (N)
    !! @param[in]  b  Right-hand side vector (N)
    !!
    !! @note A and b are not modified (copies are used internally)
    !! @note Stops execution if matrix is singular or dimensions are incompatible
    !!
    !! Algorithm:
    !! 1. Validate dimensions (A square, b and x compatible)
    !! 2. Create working copies (DGESV destroys input)
    !! 3. Call DGESV: A = P·L·U, solve L·U·x = P·b
    !! 4. Check INFO flag (0=success, >0=singular, <0=invalid argument)
    !!
    !! @see LAPACK DGESV documentation: https://netlib.org/lapack/explore-html/
    subroutine solve_linear_system(A, x, b)
        real(dp), intent(in) :: A(:,:), b(:)
        real(dp), intent(out) :: x(:)

        integer :: N, LDA, LDB, INFO, NRHS
        integer, allocatable :: IPIV(:)

        real(dp), allocatable :: A_copy(:,:), b_copy(:)

        N = size(A, 1)
        
        if (size(A, 1) /= size(A, 2)) then
            error stop "solve_linear_system: A must be square!"
        end if
        
        if (size(b) /= N .or. size(x) /= N) then
            error stop "solve_linear_system: incompatible dimensions!"
        end if

        ! LAPACK parameters
        NRHS = 1    ! Number of right-hand sides (solving A·x = b, not A·X = B)
        LDA = N     ! Leading dimension of A
        LDB = N     ! Leading dimension of b
        
        allocate(A_copy(N, N))
        allocate(b_copy(N))
        allocate(IPIV(N))

        A_copy = A
        b_copy = b

        call DGESV(N, NRHS, A_copy, LDA, IPIV, b_copy, LDB, INFO)

        if (info == 0) then
            ! Sucess
            x = b_copy
            
        else if (info > 0) then
            ! Singular matrix (there is no unique solution)
            write(*, '(A,I0,A)') "ERROR: Matrix is singular (pivot ", info, " is zero)"
            error stop "solve_linear_system: singular matrix!"
        else
            ! info < 0: invalid argument
            write(*, '(A,I0,A)') "ERROR: Invalid argument ", -info, " in DGESV"
            error stop "solve_linear_system: LAPACK error!"
        end if

        deallocate(A_copy, b_copy, IPIV)
    end subroutine

    !> Armijo backtracking line search for Newton-Raphson
    !!
    !! Finds step size α ∈ (0,1] that satisfies the Armijo condition:
    !! ||F(x + α·dx)|| < (1 - c·α)·||F(x)||
    !!
    !! This ensures sufficient decrease in the residual norm, providing
    !! global convergence guarantees for Newton's method.
    !!
    !! @param[in] x      Current iterate [k, Λ] (size N↑+M)
    !! @param[in] dx     Search direction (Newton step)
    !! @param[in] F_old  Residual at current point
    !! @param[in] I      Charge quantum numbers
    !! @param[in] J      Spin quantum numbers
    !! @param[in] L      Number of lattice sites
    !! @param[in] U      Hubbard interaction strength
    !! @return    alpha  Step size (1.0 if full step works, smaller if backtracking)
    !!
    !! Algorithm:
    !! 1. Try full step (α = 1)
    !! 2. Compute trial point: x_trial = x + α·dx
    !! 3. Evaluate residual: F_trial = F(x_trial)
    !! 4. Check Armijo condition: ||F_trial|| < (1 - c·α)·||F_old||
    !! 5. If satisfied: return α
    !! 6. If not: reduce α ← α/2 and repeat (max 20 iterations)
    !! 7. Fallback: return best α found (even if Armijo fails)
    !!
    !! @note If Armijo condition never satisfied, returns α that gives smallest residual
    !! @note Prints warning if line search fails completely
    function line_search(x, dx, F_old, I, J, L, U) result(alpha)
        real(dp), intent(in) :: x(:), dx(:), F_old(:), I(:), J(:)
        real(dp), intent(in) :: U
        integer, intent(in) :: L
        real(dp) :: alpha, best_alpha, best_norm, norm_F_old, norm_F_trial
        real(dp), allocatable :: k(:), Lambda(:), x_trial(:), F_trial(:)
        integer :: Nup, M, step

        alpha = 1.0_dp              ! Start with full Newton step
        best_alpha = alpha
        best_norm = HUGE(1.0_dp)    ! Worst possible norm (infinity)

        norm_F_old = NORM2(F_old)
        Nup = size(I)
        M = size(J)

        allocate(x_trial(Nup+M))
        allocate(k(Nup), Lambda(M))

        do step = 1, MAX_LS_ITER
            x_trial = x + alpha * dx

            k = x_trial(1:Nup)
            Lambda = x_trial(Nup+1:)

            F_trial = compute_residual(k, Lambda, I, J, L, U)
            norm_F_trial = NORM2(F_trial)

            if (norm_F_trial < best_norm) then
                best_alpha = alpha
                best_norm = norm_F_trial
            end if

            if (norm_F_trial < (1.0_dp - ARMIJO_C * alpha) * norm_F_old) then
                deallocate(x_trial, k, Lambda)
                return
            end if

            alpha = alpha / 2.0_dp
        end do

        alpha = best_alpha
        if (best_norm < norm_F_old) then
            print *, "Warning: Armijo failed, using best alpha =", alpha, &
                    " (reduction:", best_norm/norm_F_old, ")"
        else
            print *, "ERROR: Line search completely failed!"
        end if

        deallocate(x_trial, k, Lambda)
    end function

    !> Newton-Raphson solver for Bethe Ansatz equations
    !!
    !! Solves the nonlinear system F(x) = 0 where x = [k, Λ] are the
    !! charge and spin rapidities of the Lieb-Wu equations.
    !!
    !! Uses Newton-Raphson iteration with:
    !! - Analytical Jacobian (quadratic convergence)
    !! - Armijo line search (global convergence)
    !! - Special handling for U=0 (analytical solution)
    !!
    !! @param[inout] x         Solution vector [k₁,...,k_{N↑}, Λ₁,...,Λ_M]
    !!                         Input: initial guess | Output: solution (if converged)
    !! @param[in]    I         Charge quantum numbers (N↑)
    !! @param[in]    J         Spin quantum numbers (M = N↓)
    !! @param[in]    L         Number of lattice sites
    !! @param[in]    U         Hubbard interaction strength
    !! @param[out]   converged .true. if converged, .false. otherwise
    !!
    !! Algorithm:
    !! 1. Special case U≈0: return analytical solution k_j = 2π·I_j/L, Λ=0
    !! 2. Newton iteration (max 50 iterations):
    !!    a) Compute residual F(x)
    !!    b) Check convergence: ||F|| < TOL
    !!    c) Compute Jacobian J(x)
    !!    d) Solve linear system: J·dx = -F
    !!    e) Line search: find α such that ||F(x+α·dx)|| < ||F(x)||
    !!    f) Update: x ← x + α·dx
    !!    g) Check stagnation: ||dx||/||x|| < TOL
    !! 3. Return converged flag
    !!
    !! Convergence criteria:
    !! - ||F(x)|| < 1e-10 (residual tolerance)
    !! - ||dx||/||x|| < 1e-10 (relative step size)
    !!
    !! @note For U=0, returns exact solution immediately (no iteration needed)
    !! @note Prints warnings if stagnation or non-convergence detected
    !!
    !! @see compute_residual, compute_jacobian, solve_linear_system, line_search
    subroutine solve_newton(x, I, J, L, U, converged)
        real(dp), intent(inout) :: x(:)
        real(dp), intent(in) :: I(:), J(:)
        integer, intent(in) :: L
        real(dp), intent(in) :: U
        logical, intent(out) :: converged

        integer :: Nup, M, iter
        real(dp), allocatable :: k(:), Lambda(:), F(:), Jacobian(:,:), dx(:), neg_F(:)
        real(dp) :: norm_F, alpha

        Nup = size(I)
        M = size(J)

        allocate(k(Nup), Lambda(M))
        allocate(F(Nup + M), neg_F(Nup + M), dx(Nup + M))
        allocate(Jacobian(Nup + M, Nup + M))

        converged = .false.

        if (abs(U) < U_SMALL) then
            ! Analytical solution for U=0 (free Fermi gas)
            ! k_j = 2π·I_j/L
            ! Lambda is arbitrary (does not appear in the equations)
            x(1:Nup) = TWOPI * I / real(L, dp)
            x(Nup+1:) = 0.0_dp  ! Lambda = 0 (arbitrary)
            converged = .true.
            return
        end if

        do iter = 1, NEWTON_MAX_ITER
            k = x(1:Nup)
            Lambda = x(Nup+1:)

            F = compute_residual(k, Lambda, I, J, L, U)
            norm_F = NORM2(F)

            if (norm_F < NEWTON_TOL) then
                converged = .true.
                exit
            end if

            Jacobian = compute_jacobian(k, Lambda, L, U)

            neg_F = -F
            call solve_linear_system(Jacobian, dx, neg_F)

            alpha = line_search(x, dx, F, I, J, L, U)

            x = x + alpha * dx

            if (NORM2(dx) / (MAX(1.0_dp, NORM2(x))) < NEWTON_TOL) then
                if (norm_F > NEWTON_TOL) then
                    ! Stagnated without converging
                    print *, "Warning: Newton-Raphson stagnated at iteration", iter
                    print *, "  Residual norm:", norm_F
                    print *, "  Step norm:", NORM2(dx)
                else
                    ! Stagnated but converged
                    converged = .true.
                end if
                exit
            end if
        end do

        if (.not. converged) then
            print *, "Warning: Newton did not converge in", NEWTON_MAX_ITER, "iterations."
            print *, "  Final residual norm:", norm_F
        end if

        deallocate(F, Jacobian, dx, neg_F, k, Lambda)
    end subroutine

end module nonlinear_solvers