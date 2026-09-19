!> LAPACK diagonalization helpers for Kohn--Sham Hamiltonians.
module lapack_wrapper
    use lsda_constants, only: dp
    use lsda_errors, only: ERROR_INVALID_INPUT, ERROR_LAPACK_INVALID_ARG, &
                           ERROR_CONVERGENCE_FAILED, ERROR_SIZE_MISMATCH, ERROR_SUCCESS
    implicit none
    private

    type, public :: diag_workspace_t
        real(dp), allocatable :: work(:), rwork(:), d(:), e(:), real_vectors(:,:)
        complex(dp), allocatable :: cwork(:)
        integer, allocatable :: iwork(:), isuppz(:)
    end type diag_workspace_t

    public :: validate_diagonalization_inputs
    public :: diagonalize_symmetric_real, diagonalize_symmetric_real_values_only
    public :: diagonalize_hermitian_complex, diagonalize_hermitian_complex_values_only
    public :: diagonalize_open_tridiagonal, diagonalize_open_tridiagonal_complex, diagonalize_symmetric_real_partial
    public :: diagonalize_hermitian_complex_partial, cleanup_diag_workspace

    interface
        subroutine DSTEVR(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info)
            import :: dp
            character(len=1), intent(in) :: jobz, range
            integer, intent(in) :: n, il, iu, ldz, lwork, liwork
            real(dp), intent(inout) :: d(*), e(*)
            real(dp), intent(in) :: vl, vu, abstol
            integer, intent(out) :: m, isuppz(*)
            real(dp), intent(out) :: w(*), z(ldz,*), work(*)
            integer, intent(out) :: iwork(*), info
        end subroutine DSTEVR
        subroutine DSYEVR(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info)
            import :: dp
            character(len=1), intent(in) :: jobz, range, uplo
            integer, intent(in) :: n, lda, il, iu, ldz, lwork, liwork
            real(dp), intent(inout) :: a(lda,*)
            real(dp), intent(in) :: vl, vu, abstol
            integer, intent(out) :: m, isuppz(*)
            real(dp), intent(out) :: w(*), z(ldz,*), work(*)
            integer, intent(out) :: iwork(*), info
        end subroutine DSYEVR
        subroutine ZHEEVR(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)
            import :: dp
            character(len=1), intent(in) :: jobz, range, uplo
            integer, intent(in) :: n, lda, il, iu, ldz, lwork, lrwork, liwork
            complex(dp), intent(inout) :: a(lda,*), z(ldz,*), work(*)
            real(dp), intent(in) :: vl, vu, abstol
            integer, intent(out) :: m, isuppz(*)
            real(dp), intent(out) :: w(*), rwork(*)
            integer, intent(out) :: iwork(*), info
        end subroutine ZHEEVR
    end interface

contains

    !> Validate legacy full-spectrum real diagonalization arrays.
    subroutine validate_diagonalization_inputs(L, H, eigvals, eigvecs, compute_vectors, ierr)
        integer, intent(in) :: L
        real(dp), intent(in) :: H(:,:), eigvals(:)
        real(dp), intent(in), optional :: eigvecs(:,:)
        logical, intent(in) :: compute_vectors
        integer, intent(out) :: ierr

        ierr = ERROR_SUCCESS
        if (L <= 0) then; ierr = ERROR_INVALID_INPUT; return; end if
        if (size(H, 1) /= L .or. size(H, 2) /= L .or. size(eigvals) /= L) then
            ierr = ERROR_SIZE_MISMATCH; return
        end if
        if (compute_vectors) then
            if (.not. present(eigvecs)) then; ierr = ERROR_INVALID_INPUT; return; end if
            if (size(eigvecs, 1) /= L .or. size(eigvecs, 2) /= L) ierr = ERROR_SIZE_MISMATCH
        end if
    end subroutine validate_diagonalization_inputs

    !> Release reusable diagonalization workspace.
    subroutine cleanup_diag_workspace(workspace)
        type(diag_workspace_t), intent(inout) :: workspace
        if (allocated(workspace%work)) deallocate(workspace%work)
        if (allocated(workspace%rwork)) deallocate(workspace%rwork)
        if (allocated(workspace%cwork)) deallocate(workspace%cwork)
        if (allocated(workspace%iwork)) deallocate(workspace%iwork)
        if (allocated(workspace%isuppz)) deallocate(workspace%isuppz)
        if (allocated(workspace%d)) deallocate(workspace%d)
        if (allocated(workspace%e)) deallocate(workspace%e)
        if (allocated(workspace%real_vectors)) deallocate(workspace%real_vectors)
    end subroutine cleanup_diag_workspace

    !> Ensure that a real MRRR workspace can solve a system of size L.
    subroutine ensure_real_workspace(workspace, L)
        type(diag_workspace_t), intent(inout) :: workspace
        integer, intent(in) :: L
        if (.not. allocated(workspace%work) .or. size(workspace%work) /= 26 * L) then
            if (allocated(workspace%work)) deallocate(workspace%work)
            allocate(workspace%work(max(1, 26 * L)))
        end if
        if (.not. allocated(workspace%iwork) .or. size(workspace%iwork) /= 10 * L) then
            if (allocated(workspace%iwork)) deallocate(workspace%iwork)
            allocate(workspace%iwork(max(1, 10 * L)))
        end if
        if (.not. allocated(workspace%isuppz) .or. size(workspace%isuppz) /= 2 * L) then
            if (allocated(workspace%isuppz)) deallocate(workspace%isuppz)
            allocate(workspace%isuppz(max(2, 2 * L)))
        end if
        if (.not. allocated(workspace%d) .or. size(workspace%d) /= L) then
            if (allocated(workspace%d)) deallocate(workspace%d)
            allocate(workspace%d(L))
        end if
        if (.not. allocated(workspace%e) .or. size(workspace%e) /= max(1, L - 1)) then
            if (allocated(workspace%e)) deallocate(workspace%e)
            allocate(workspace%e(max(1, L - 1)))
        end if
    end subroutine ensure_real_workspace

    !> Ensure that a complex Hermitian MRRR workspace can solve size L.
    subroutine ensure_complex_workspace(workspace, L)
        type(diag_workspace_t), intent(inout) :: workspace
        integer, intent(in) :: L
        call ensure_real_workspace(workspace, L)
        if (.not. allocated(workspace%cwork) .or. size(workspace%cwork) /= 2 * L) then
            if (allocated(workspace%cwork)) deallocate(workspace%cwork)
            allocate(workspace%cwork(max(1, 2 * L)))
        end if
        if (.not. allocated(workspace%rwork) .or. size(workspace%rwork) /= 24 * L) then
            if (allocated(workspace%rwork)) deallocate(workspace%rwork)
            allocate(workspace%rwork(max(1, 24 * L)))
        end if
    end subroutine ensure_complex_workspace

    !> Diagonalize an open-chain tridiagonal Hamiltonian using DSTEVR.
    !! @param[in] potential On-site effective potential (length L)
    !! @param[in] n_vec Number of lowest eigenpairs requested
    !! @param[out] eigvals Lowest eigenvalues (length >= n_vec)
    !! @param[out] eigvecs Corresponding eigenvectors (L x n_vec)
    !! @param[inout] workspace Reused LAPACK scratch storage
    !! @param[out] ierr Error code
    subroutine diagonalize_open_tridiagonal(potential, L, n_vec, eigvals, eigvecs, workspace, ierr)
        real(dp), intent(in) :: potential(:)
        integer, intent(in) :: L, n_vec
        real(dp), intent(out) :: eigvals(:), eigvecs(:,:)
        type(diag_workspace_t), intent(inout) :: workspace
        integer, intent(out) :: ierr
        integer :: info, m
        ierr = ERROR_SUCCESS
        if (L <= 0 .or. n_vec < 1 .or. n_vec > L) then; ierr = ERROR_INVALID_INPUT; return; end if
        if (size(potential) /= L .or. size(eigvals) < n_vec .or. size(eigvecs,1) /= L .or. size(eigvecs,2) < n_vec) then
            ierr = ERROR_SIZE_MISMATCH; return
        end if
        call ensure_real_workspace(workspace, L)
        workspace%d(1:L) = potential
        if (L > 1) workspace%e(1:L-1) = -1.0_dp
        call DSTEVR('V', 'I', L, workspace%d, workspace%e, 0.0_dp, 0.0_dp, 1, n_vec, 0.0_dp, m, eigvals, eigvecs, L, &
                    workspace%isuppz, workspace%work, size(workspace%work), workspace%iwork, size(workspace%iwork), info)
        if (info < 0 .or. m /= n_vec) then
            ierr = ERROR_LAPACK_INVALID_ARG
        else if (info > 0) then
            ierr = ERROR_CONVERGENCE_FAILED
        end if
    end subroutine diagonalize_open_tridiagonal

    !> Diagonalize an open-chain Hamiltonian and return complex eigenvectors.
    !!
    !! Open-boundary hopping and on-site potentials are real, so DSTEVR is
    !! still used; its real eigenvectors are promoted without changing phase.
    subroutine diagonalize_open_tridiagonal_complex(potential, L, n_vec, eigvals, eigvecs, workspace, ierr)
        real(dp), intent(in) :: potential(:)
        integer, intent(in) :: L, n_vec
        real(dp), intent(out) :: eigvals(:)
        complex(dp), intent(out) :: eigvecs(:,:)
        type(diag_workspace_t), intent(inout) :: workspace
        integer, intent(out) :: ierr

        real(dp), allocatable :: real_vectors(:,:)

        if (.not. allocated(workspace%real_vectors) .or. size(workspace%real_vectors, 1) /= L .or. &
            size(workspace%real_vectors, 2) /= n_vec) then
            if (allocated(workspace%real_vectors)) deallocate(workspace%real_vectors)
            allocate(workspace%real_vectors(L, n_vec))
        end if

        ! Detach the buffer for the duration of the call. Passing
        ! `workspace%real_vectors` as the definable `eigvecs` dummy while
        ! `workspace` itself is associated with the definable `workspace` dummy
        ! places the same storage under two definable dummies, one of them a
        ! subobject of the other's actual argument, which F2018 15.5.2.13
        ! forbids. It happens to work today only because `ensure_real_workspace`
        ! never touches `real_vectors`; `move_alloc` removes the aliasing rather
        ! than relying on that. Both moves are O(1), no array is copied.
        !
        ! If the callee is ever changed to allocate `workspace%real_vectors`
        ! itself, the second `move_alloc` would discard that allocation; use a
        ! separate local buffer and copy back instead.
        call move_alloc(workspace%real_vectors, real_vectors)
        call diagonalize_open_tridiagonal(potential, L, n_vec, eigvals, real_vectors, workspace, ierr)
        if (ierr == ERROR_SUCCESS) eigvecs(:,1:n_vec) = cmplx(real_vectors(:,1:n_vec), 0.0_dp, kind=dp)
        call move_alloc(real_vectors, workspace%real_vectors)
    end subroutine diagonalize_open_tridiagonal_complex

    !> Diagonalize only the lowest real-symmetric eigenpairs with DSYEVR.
    subroutine diagonalize_symmetric_real_partial(H, L, n_vec, eigvals, eigvecs, workspace, ierr)
        real(dp), intent(inout) :: H(:,:)
        integer, intent(in) :: L, n_vec
        real(dp), intent(out) :: eigvals(:), eigvecs(:,:)
        type(diag_workspace_t), intent(inout) :: workspace
        integer, intent(out) :: ierr
        integer :: info, m
        ierr = ERROR_SUCCESS
        if (L <= 0 .or. n_vec < 1 .or. n_vec > L) then; ierr = ERROR_INVALID_INPUT; return; end if
        if (size(H,1) /= L .or. size(H,2) /= L .or. size(eigvals) < n_vec .or. size(eigvecs,1) /= L .or. size(eigvecs,2) < n_vec) then
            ierr = ERROR_SIZE_MISMATCH; return
        end if
        call ensure_real_workspace(workspace, L)
        call DSYEVR('V', 'I', 'U', L, H, L, 0.0_dp, 0.0_dp, 1, n_vec, 0.0_dp, m, eigvals, eigvecs, L, workspace%isuppz, &
                    workspace%work, size(workspace%work), workspace%iwork, size(workspace%iwork), info)
        if (info < 0 .or. m /= n_vec) then; ierr = ERROR_LAPACK_INVALID_ARG
        else if (info > 0) then; ierr = ERROR_CONVERGENCE_FAILED
        end if
    end subroutine diagonalize_symmetric_real_partial

    !> Diagonalize only the lowest complex-Hermitian eigenpairs with ZHEEVR.
    subroutine diagonalize_hermitian_complex_partial(H, L, n_vec, eigvals, eigvecs, workspace, ierr)
        complex(dp), intent(inout) :: H(:,:)
        integer, intent(in) :: L, n_vec
        real(dp), intent(out) :: eigvals(:)
        complex(dp), intent(out) :: eigvecs(:,:)
        type(diag_workspace_t), intent(inout) :: workspace
        integer, intent(out) :: ierr
        integer :: info, m
        ierr = ERROR_SUCCESS
        if (L <= 0 .or. n_vec < 1 .or. n_vec > L) then; ierr = ERROR_INVALID_INPUT; return; end if
        if (size(H,1) /= L .or. size(H,2) /= L .or. size(eigvals) < n_vec .or. size(eigvecs,1) /= L .or. size(eigvecs,2) < n_vec) then
            ierr = ERROR_SIZE_MISMATCH; return
        end if
        call ensure_complex_workspace(workspace, L)
        call ZHEEVR('V', 'I', 'U', L, H, L, 0.0_dp, 0.0_dp, 1, n_vec, 0.0_dp, m, eigvals, eigvecs, L, workspace%isuppz, &
                    workspace%cwork, size(workspace%cwork), workspace%rwork, size(workspace%rwork), workspace%iwork, size(workspace%iwork), info)
        if (info < 0 .or. m /= n_vec) then; ierr = ERROR_LAPACK_INVALID_ARG
        else if (info > 0) then; ierr = ERROR_CONVERGENCE_FAILED
        end if
    end subroutine diagonalize_hermitian_complex_partial

    !> Backward-compatible full real diagonalization.
    subroutine diagonalize_symmetric_real(H, L, eigvals, eigvecs, ierr)
        real(dp), intent(in) :: H(:,:)
        integer, intent(in) :: L
        real(dp), intent(out) :: eigvals(:), eigvecs(:,:)
        integer, intent(out) :: ierr
        real(dp), allocatable :: H_work(:,:)
        type(diag_workspace_t) :: workspace
        call validate_diagonalization_inputs(L, H, eigvals, eigvecs, .true., ierr)
        if (ierr /= ERROR_SUCCESS) return
        allocate(H_work(L,L)); H_work = H
        call diagonalize_symmetric_real_partial(H_work, L, L, eigvals, eigvecs, workspace, ierr)
        call cleanup_diag_workspace(workspace); deallocate(H_work)
    end subroutine diagonalize_symmetric_real

    !> Backward-compatible real eigenvalue-only diagonalization.
    subroutine diagonalize_symmetric_real_values_only(H, L, eigvals, ierr)
        real(dp), intent(in) :: H(:,:)
        integer, intent(in) :: L
        real(dp), intent(out) :: eigvals(:)
        integer, intent(out) :: ierr
        real(dp), allocatable :: vectors(:,:), H_work(:,:)
        type(diag_workspace_t) :: workspace
        call validate_diagonalization_inputs(L, H, eigvals, compute_vectors=.false., ierr=ierr)
        if (ierr /= ERROR_SUCCESS) return
        allocate(H_work(L,L), vectors(L,L)); H_work = H
        call diagonalize_symmetric_real_partial(H_work, L, L, eigvals, vectors, workspace, ierr)
        call cleanup_diag_workspace(workspace); deallocate(H_work, vectors)
    end subroutine diagonalize_symmetric_real_values_only

    !> Backward-compatible full complex-Hermitian diagonalization.
    subroutine diagonalize_hermitian_complex(H, L, eigvals, eigvecs, ierr)
        complex(dp), intent(in) :: H(:,:)
        integer, intent(in) :: L
        real(dp), intent(out) :: eigvals(:)
        complex(dp), intent(out) :: eigvecs(:,:)
        integer, intent(out) :: ierr
        complex(dp), allocatable :: H_work(:,:)
        type(diag_workspace_t) :: workspace
        ierr = ERROR_SUCCESS
        if (L <= 0) then; ierr = ERROR_INVALID_INPUT; return; end if
        if (size(H,1) /= L .or. size(H,2) /= L .or. size(eigvals) /= L .or. size(eigvecs,1) /= L .or. size(eigvecs,2) /= L) then; ierr = ERROR_SIZE_MISMATCH; return; end if
        allocate(H_work(L,L)); H_work = H
        call diagonalize_hermitian_complex_partial(H_work, L, L, eigvals, eigvecs, workspace, ierr)
        call cleanup_diag_workspace(workspace); deallocate(H_work)
    end subroutine diagonalize_hermitian_complex

    !> Backward-compatible complex-Hermitian eigenvalue-only diagonalization.
    subroutine diagonalize_hermitian_complex_values_only(H, L, eigvals, ierr)
        complex(dp), intent(in) :: H(:,:)
        integer, intent(in) :: L
        real(dp), intent(out) :: eigvals(:)
        integer, intent(out) :: ierr
        complex(dp), allocatable :: H_work(:,:), vectors(:,:)
        type(diag_workspace_t) :: workspace
        ierr = ERROR_SUCCESS
        if (L <= 0) then; ierr = ERROR_INVALID_INPUT; return; end if
        if (size(H,1) /= L .or. size(H,2) /= L .or. size(eigvals) /= L) then; ierr = ERROR_SIZE_MISMATCH; return; end if
        allocate(H_work(L,L), vectors(L,L)); H_work = H
        call diagonalize_hermitian_complex_partial(H_work, L, L, eigvals, vectors, workspace, ierr)
        call cleanup_diag_workspace(workspace); deallocate(H_work, vectors)
    end subroutine diagonalize_hermitian_complex_values_only
end module lapack_wrapper
