module nonlinear_solvers
    use lsda_constants, only: dp
    use bethe_equations
    implicit none
    private

    public :: solve_newton
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

end module nonlinear_solvers