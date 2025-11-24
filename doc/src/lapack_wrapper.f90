module lapack_wrapper
    use lsda_constants, only: dp
    use lsda_errors, only: ERROR_INVALID_INPUT, ERROR_LAPACK_FAILED, &
                        ERROR_SIZE_MISMATCH, ERROR_LAPACK_INVALID_ARG, &
                        ERROR_CONVERGENCE_FAILED, ERROR_SUCCESS
    implicit none
    private

    public :: validate_diagonalization_inputs
    public :: diagonalize_symmetric_real
    public :: diagonalize_symmetric_real_values_only
    public :: diagonalize_hermitian_complex
    public :: diagonalize_hermitian_complex_values_only

    interface
        subroutine DSYEVD(jobz, uplo, n, a, lda, w, &
                        work, lwork, iwork, liwork, info)
            import :: dp
            character(len=1), intent(in) :: jobz
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n
            integer, intent(in) :: lda
            real(dp), intent(inout) :: a(lda,*)
            real(dp), intent(out) :: w(*)
            real(dp), intent(out) :: work(*)
            integer, intent(in) :: lwork
            integer, intent(out) :: iwork(*)
            integer, intent(in) :: liwork
            integer, intent(out) :: info
        end subroutine DSYEVD

        subroutine ZHEEVD(jobz, uplo, n, a, lda, w, &
                        work, lwork, rwork, lrwork, &
                        iwork, liwork, info)
            import :: dp
            character(len=1), intent(in) :: jobz
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n
            integer, intent(in) :: lda
            complex(dp), intent(inout) :: a(lda,*)
            real(dp), intent(out) :: w(*)
            complex(dp), intent(out) :: work(*)
            integer, intent(in) :: lwork
            real(dp), intent(out) :: rwork(*)
            integer, intent(in) :: lrwork
            integer, intent(out) :: iwork(*)
            integer, intent(in) :: liwork
            integer, intent(out) :: info
        end subroutine ZHEEVD
    end interface

contains

    subroutine validate_diagonalization_inputs(L, H, eigvals, eigvecs, compute_vectors, ierr)
        integer, intent(in) :: L
        real(dp), intent(in) :: H(:,:), eigvals(:)
        real(dp), intent(in), optional :: eigvecs(:,:)
        logical, intent(in) :: compute_vectors
        integer, intent(out) :: ierr

        ierr = 0

        if (L <= 0) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (size(H, 1) /= L .or. size(H, 2) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        if (size(eigvals) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        if (compute_vectors) then
            if (.not. present(eigvecs)) then
                ierr = ERROR_INVALID_INPUT
                return
            end if

            if (size(eigvecs, 1) /= L .or. size(eigvecs, 2) /= L) then
                ierr = ERROR_SIZE_MISMATCH
                return
            end if
        end if
    end subroutine validate_diagonalization_inputs

    !> @brief Diagonalize real symmetric matrix (eigenvalues + eigenvectors)
    !!
    !! Solves H*v = lambda*v for symmetric H using LAPACK DSYEVD.
    !! Returns eigenvalues in ascending order and eigenvectors as columns.
    !!
    !! @param[in] H Real symmetric matrix (L x L)
    !! @param[in] L Dimension of matrix
    !! @param[out] eigvals Eigenvalues in ascending order (length L)
    !! @param[out] eigvecs Eigenvectors as columns (L x L)
    !! @param[out] ierr Error code (0 = success)
    subroutine diagonalize_symmetric_real(H, L, eigvals, eigvecs, ierr)
        real(dp), intent(in) :: H(:,:)
        integer, intent(in) :: L
        real(dp), intent(out) :: eigvals(:)
        real(dp), intent(out) :: eigvecs(:,:)
        integer, intent(out) :: ierr
        
        real(dp), allocatable :: H_work(:,:)
        real(dp), allocatable :: work(:)
        integer, allocatable :: iwork(:)
        integer :: lwork, liwork, info

        call validate_diagonalization_inputs(L, H, eigvals, eigvecs, .true., ierr)

        if (ierr /= ERROR_SUCCESS) return

        allocate(H_work(L,L))
        H_work = H

        lwork = -1
        liwork = -1
        allocate(work(1), iwork(1))

        call DSYEVD('V', 'U', L, H_work, L, eigvals, work, &
                lwork, iwork, liwork, info)

        if (info /= 0) then
            ierr = ERROR_LAPACK_INVALID_ARG
            deallocate(H_work, work, iwork)
            return
        end if

        lwork = int(work(1))
        liwork = iwork(1)
        deallocate(work, iwork)
        allocate(work(lwork), iwork(liwork))

        ! Actual diagonalization
        call DSYEVD('V', 'U', L, H_work, L, eigvals, work, &
                lwork, iwork, liwork, info)

        if (info < 0) then
            ierr = ERROR_LAPACK_INVALID_ARG
        else if (info > 0) then
            ierr = ERROR_CONVERGENCE_FAILED
        else
            ierr = ERROR_SUCCESS
            eigvecs = H_work
        end if

        deallocate(H_work, work, iwork)
    end subroutine diagonalize_symmetric_real

    !> @brief Diagonalize real symmetric matrix (eigenvalues only)
    !!
    !! Faster version that only computes eigenvalues, not eigenvectors.
    !!
    !! @param[in] H Real symmetric matrix (L x L)
    !! @param[in] L Dimension of matrix
    !! @param[out] eigvals Eigenvalues in ascending order (length L)
    !! @param[out] ierr Error code (0 = success)
    subroutine diagonalize_symmetric_real_values_only(H, L, eigvals, ierr)
        real(dp), intent(in) :: H(:,:)
        integer, intent(in) :: L
        real(dp), intent(out) :: eigvals(:)
        integer, intent(out) :: ierr
        
        real(dp), allocatable :: H_work(:,:)
        real(dp), allocatable :: work(:)
        integer, allocatable :: iwork(:)
        integer :: lwork, liwork, info

        call validate_diagonalization_inputs(L, H, eigvals, compute_vectors=.false., ierr=ierr)

        if (ierr /= ERROR_SUCCESS) return

        allocate(H_work(L,L))
        H_work = H

        lwork = -1
        liwork = -1
        allocate(work(1), iwork(1))

        call DSYEVD('N', 'U', L, H_work, L, eigvals, work, &
                lwork, iwork, liwork, info)

        if (info /= 0) then
            ierr = ERROR_LAPACK_INVALID_ARG
            deallocate(H_work, work, iwork)
            return
        end if

        lwork = int(work(1))
        liwork = iwork(1)
        deallocate(work, iwork)
        allocate(work(lwork), iwork(liwork))

        ! Actual diagonalization
        call DSYEVD('N', 'U', L, H_work, L, eigvals, work, &
                    lwork, iwork, liwork, info)

        if (info < 0) then
            ierr = ERROR_LAPACK_INVALID_ARG
        else if (info > 0) then
            ierr = ERROR_CONVERGENCE_FAILED
        else
            ierr = ERROR_SUCCESS
        end if

        deallocate(H_work, work, iwork)
    end subroutine diagonalize_symmetric_real_values_only

    !> @brief Diagonalize complex Hermitian matrix (eigenvalues + eigenvectors)
    !!
    !! Solves H*v = lambda*v for Hermitian H using LAPACK ZHEEVD.
    !! Eigenvalues are real even though H is complex.
    !!
    !! @param[in] H Complex Hermitian matrix (L x L)
    !! @param[in] L Dimension of matrix
    !! @param[out] eigvals Real eigenvalues in ascending order (length L)
    !! @param[out] eigvecs Complex eigenvectors as columns (L x L)
    !! @param[out] ierr Error code (0 = success)
    subroutine diagonalize_hermitian_complex(H, L, eigvals, eigvecs, ierr)
        complex(dp), intent(in) :: H(:,:)
        integer, intent(in) :: L
        real(dp), intent(out) :: eigvals(:)
        complex(dp), intent(out) :: eigvecs(:,:)
        integer, intent(out) :: ierr
        
        complex(dp), allocatable :: H_work(:,:)
        complex(dp), allocatable :: work(:)
        real(dp), allocatable :: rwork(:)
        integer, allocatable :: iwork(:)
        integer :: lwork, lrwork, liwork, info

        ! Validate inputs
        ierr = ERROR_SUCCESS
        if (L <= 0) then
            ierr = ERROR_INVALID_INPUT
            return
        end if
        if (size(H, 1) /= L .or. size(H, 2) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if
        if (size(eigvals) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if
        if (size(eigvecs, 1) /= L .or. size(eigvecs, 2) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        allocate(H_work(L,L))
        H_work = H

        lwork = -1
        lrwork = -1
        liwork = -1
        allocate(work(1), rwork(1), iwork(1))

        call ZHEEVD('V', 'U', L, H_work, L, eigvals, work, &
                lwork, rwork, lrwork, iwork, liwork, info)

        if (info /= 0) then
            ierr = ERROR_LAPACK_INVALID_ARG
            deallocate(H_work, work, rwork, iwork)
            return
        end if

        lwork = int(real(work(1), kind=dp))
        lrwork = int(rwork(1))
        liwork = iwork(1)
        deallocate(work, rwork, iwork)
        allocate(work(lwork), iwork(liwork), rwork(lrwork))

        ! Actual diagonalization
        call ZHEEVD('V', 'U', L, H_work, L, eigvals, work, &
                lwork, rwork, lrwork, iwork, liwork, info)

        if (info < 0) then
            ierr = ERROR_LAPACK_INVALID_ARG
        else if (info > 0) then
            ierr = ERROR_CONVERGENCE_FAILED
        else
            ierr = ERROR_SUCCESS
            eigvecs = H_work
        end if

        deallocate(H_work, work, rwork, iwork)
    end subroutine diagonalize_hermitian_complex

    subroutine diagonalize_hermitian_complex_values_only(H, L, eigvals, ierr)
        complex(dp), intent(in) :: H(:,:)
        integer, intent(in) :: L
        real(dp), intent(out) :: eigvals(:)
        integer, intent(out) :: ierr
        
        complex(dp), allocatable :: H_work(:,:)
        complex(dp), allocatable :: work(:)
        real(dp), allocatable :: rwork(:)
        integer, allocatable :: iwork(:)
        integer :: lwork, lrwork, liwork, info

        ! Validate inputs
        ierr = ERROR_SUCCESS
        if (L <= 0) then
            ierr = ERROR_INVALID_INPUT
            return
        end if
        if (size(H, 1) /= L .or. size(H, 2) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if
        if (size(eigvals) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        allocate(H_work(L,L))
        H_work = H

        lwork = -1
        lrwork = -1
        liwork = -1
        allocate(work(1), rwork(1), iwork(1))

        call ZHEEVD('N', 'U', L, H_work, L, eigvals, work, &
                lwork, rwork, lrwork, iwork, liwork, info)

        if (info /= 0) then
            ierr = ERROR_LAPACK_INVALID_ARG
            deallocate(H_work, work, rwork, iwork)
            return
        end if

        lwork = int(real(work(1), kind=dp))
        lrwork = int(rwork(1))
        liwork = iwork(1)
        deallocate(work, rwork, iwork)
        allocate(work(lwork), iwork(liwork), rwork(lrwork))

        ! Actual diagonalization
        call ZHEEVD('N', 'U', L, H_work, L, eigvals, work, &
                lwork, rwork, lrwork, iwork, liwork, info)

        if (info < 0) then
            ierr = ERROR_LAPACK_INVALID_ARG
        else if (info > 0) then
            ierr = ERROR_CONVERGENCE_FAILED
        else
            ierr = ERROR_SUCCESS
        end if

        deallocate(H_work, work, rwork, iwork)
    end subroutine diagonalize_hermitian_complex_values_only
end module lapack_wrapper