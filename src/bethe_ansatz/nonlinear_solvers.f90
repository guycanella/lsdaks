module nonlinear_solvers
    use lsda_constants, only: dp
    use bethe_equations
    implicit none
    private

    ! public :: solve_newton
    public :: solve_linear_system

    real(dp), parameter :: ARMIJO_C = 1.0e-4_dp
    integer, parameter :: MAX_LS_ITER = 20

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

    subroutine solve_linear_system(A, b, x)
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

        NRHS = 1
        LDA = N
        LDB = N
        
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

    function line_search(x, dx, F_old, I, J, L, U) result(alpha)
        real(dp) :: alpha, norm_F_old, norm_F_trial, x_trial(size(x)), F_trial(:)
        real(dp), allocatable :: k(:), Lambda(:)
        integer :: Nup, M, step
        real(dp), intent(in) :: x(:), dx(:), F_old(:), I(:), J(:)
        real(dp), intent(in) :: U
        integer, intent(in) :: L

        alpha = 1.0_dp
        norm_F_old = NORM2(F_old)
        Nup = size(I)
        M = size(J)

        allocate(k(Nup), Lambda(M))

        do step = 1, MAX_LS_ITER
            x_trial = x + alpha * dx

            k = x_trial(1:Nup)
            Lambda = x_trial(Nup+1:)
            
            F_trial = compute_residual(k, Lambda, I, J, L, U)
            norm_F_trial = NORM2(F_trial)

            if (norm_F_trial < (1.0_dp - ARMIJO_C * alpha) * norm_F_old) then
                deallocate(k, Lambda)
                return
            end if

            alpha = alpha / 2.0_dp
        end do

        deallocate(k, Lambda)
    end function

end module nonlinear_solvers