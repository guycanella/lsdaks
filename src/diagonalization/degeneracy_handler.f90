module degeneracy_handler
    use lsda_constants, only: dp
    use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT, &
                        ERROR_SIZE_MISMATCH, ERROR_LAPACK_FAILED, &
                        ERROR_LINEAR_DEPENDENCE
    implicit none
    private

    real(dp), parameter :: DEG_TOL = 1.0e-8_dp

    interface
        subroutine DGEQRF(M, N, A, LDA, TAU, WORK, LWORK, INFO)
            import :: dp
            integer, intent(in) :: M, N, LDA, LWORK
            real(dp), intent(inout) :: A(LDA, *)
            real(dp), intent(out) :: TAU(*)
            real(dp), intent(inout) :: WORK(*)
            integer, intent(out) :: INFO
        end subroutine DGEQRF

        subroutine DORGQR(M, N, K, A, LDA, TAU, WORK, LWORK, INFO)
            import :: dp
            integer, intent(in) :: M, N, K, LDA, LWORK
            real(dp), intent(inout) :: A(LDA, *)
            real(dp), intent(in) :: TAU(*)
            real(dp), intent(inout) :: WORK(*)
            integer, intent(out) :: INFO
        end subroutine DORGQR

        subroutine DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, &
                                                LDB, BETA, C, LDC)
            import :: dp
            character(len=1), intent(in) :: TRANSA, TRANSB
            integer, intent(in) :: M, N, K, LDA, LDB, LDC
            real(dp), intent(in) :: ALPHA, BETA
            real(dp), intent(in) :: A(LDA, *), B(LDB, *)
            real(dp), intent(inout) :: C(LDC, *)
        end subroutine DGEMM
    end interface

    public :: find_degenerate_subspaces
    public :: orthonormalize_degenerate_subspace
    public :: orthonormalize_degenerate_subspace_complex
    public :: compute_degeneracy_count
    public :: verify_orthonormality

contains

    !> @brief Find all degenerate subspaces in eigenvalue spectrum
    !!
    !! Scans eigenvalues (assumed sorted) and identifies groups where
    !! |λ_i - λ_j| < tol. Returns indices of degenerate subspaces.
    !!
    !! @param[in] eigvals Eigenvalues (sorted in ascending order)
    !! @param[in] L Number of eigenvalues
    !! @param[out] subspaces Indices of degenerate eigenvalues (n_subspaces × max_deg)
    !! @param[out] subspace_sizes Size of each degenerate subspace
    !! @param[out] n_subspaces Number of degenerate subspaces found
    !! @param[out] ierr Error code
    subroutine find_degenerate_subspaces(eigvals, L, subspaces, subspaces_sizes, &
                                                                n_subspaces, ierr)
        real(dp), intent(in) ::eigvals(:)
        integer, intent(in) :: L
        integer, allocatable, intent(out) :: subspaces(:,:)
        integer, allocatable, intent(out) :: subspaces_sizes(:)
        integer, intent(out) :: ierr
        integer, intent(out) :: n_subspaces

        integer :: i, j, deg_count, idx, max_deg

        ierr = ERROR_SUCCESS

        if (L <= 0) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (size(eigvals) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        i = 1
        n_subspaces = 0
        max_deg = 1

        do while (i <= L)
            deg_count = 1
            do j = i + 1, L
                if (abs(eigvals(j) - eigvals(i)) < DEG_TOL) then
                    deg_count = deg_count + 1
                else
                    exit
                end if
            end do

            if (deg_count > 1) then
                n_subspaces = n_subspaces + 1
                max_deg = max(max_deg, deg_count)
            end if

            i = i + deg_count
        end do

        ! If no degeneracies found, return empty arrays
        if (n_subspaces == 0) then
            allocate(subspaces(0, 0))
            allocate(subspaces_sizes(0))
            return
        end if

        allocate(subspaces(n_subspaces, max_deg))
        allocate(subspaces_sizes(n_subspaces))
        subspaces = 0

        idx = 0
        i = 1
        do while (i <= L)
            deg_count = 1
            do j = i + 1, L
                if (abs(eigvals(j) - eigvals(i)) < DEG_TOL) then
                    deg_count = deg_count + 1
                else
                    exit
                end if
            end do

            if (deg_count > 1) then
                idx = idx + 1
                subspaces_sizes(idx) = deg_count

                do j = 1, deg_count
                    subspaces(idx, j) = i + j - 1
                end do
            end if

            i = i + deg_count
        end do
    end subroutine find_degenerate_subspaces

    !> @brief Re-orthonormalize degenerate subspace using LAPACK QR
    !!
    !! Uses QR decomposition to obtain orthonormal basis for a degenerate
    !! subspace. More numerically stable than Gram-Schmidt.
    !!
    !! @param[inout] eigvecs Eigenvector matrix (columns are eigenvectors)
    !! @param[in] L Dimension of eigenvectors
    !! @param[in] indices Indices of degenerate eigenvectors to orthonormalize
    !! @param[in] n_deg Number of degenerate eigenvectors
    !! @param[out] ierr Error code
    subroutine orthonormalize_degenerate_subspace(eigvecs, L, indexes, n_deg, ierr)
        real(dp), intent(inout) :: eigvecs(:,:)
        integer, intent(in) :: L, indexes(:), n_deg
        integer, intent(out) :: ierr

        real(dp), allocatable :: A(:,:), tau(:), work(:)
        integer :: lwork, info, i

        if (n_deg <= 1) then
            ierr = ERROR_SUCCESS ! Nothing to do for a single vector
            return
        end if

        if (L <= 0) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (size(indexes) < n_deg) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        if (size(eigvecs, 1) /= L .or. size(eigvecs, 2) < L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        do i = 1, n_deg
            if (indexes(i) < 1 .or. indexes(i) > L) then
                ierr = ERROR_INVALID_INPUT
                return
            end if
        end do

        allocate(A(L, n_deg))
        do i = 1, n_deg
            A(:, i) = eigvecs(:, indexes(i))
        end do

        allocate(tau(n_deg))

        lwork = -1
        allocate(work(1))
        call DGEQRF(L, n_deg, A, L, tau, work, lwork, info)
        if (info /= 0) then
            ierr = ERROR_LAPACK_FAILED
            deallocate(A, tau, work)
            return
        end if

        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        call DGEQRF(L, n_deg, A, L, tau, work, lwork, info)
        if (info /= 0) then
            ierr = ERROR_LAPACK_FAILED
            deallocate(A, tau, work)
            return
        end if

        ! Workspace query for DORGQR
        deallocate(work)
        allocate(work(1))
        lwork = -1
        call DORGQR(L, n_deg, n_deg, A, L, tau, work, lwork, info)
        if (info /= 0) then
            ierr = ERROR_LAPACK_FAILED
            deallocate(A, tau, work)
            return
        end if

        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        call DORGQR(L, n_deg, n_deg, A, L, tau, work, lwork, info)
        if (info /= 0) then
            ierr = ERROR_LAPACK_FAILED
            deallocate(A, tau, work)
            return
        end if
        
        do i = 1, n_deg
            eigvecs(:, indexes(i)) = A(:, i)
        end do
        
        deallocate(A, tau, work)
    end subroutine orthonormalize_degenerate_subspace

    !> @brief Re-orthonormalize complex degenerate subspace using Gram-Schmidt
    !!
    !! Modified Gram-Schmidt for complex vectors. Less stable than QR
    !! but simpler for complex case.
    !!
    !! @param[inout] eigvecs Complex eigenvector matrix
    !! @param[in] L Dimension of eigenvectors
    !! @param[in] indexes Indices of degenerate eigenvectors
    !! @param[in] n_deg Number of degenerate eigenvectors
    !! @param[out] ierr Error code
    subroutine orthonormalize_degenerate_subspace_complex(eigvecs, L, indexes, n_deg, ierr)
        complex(dp), intent(inout) :: eigvecs(:,:)
        integer, intent(in) :: L, indexes(:), n_deg
        integer, intent(out) :: ierr

        integer :: j, k, idx_j, idx_k
        complex(dp) :: proj
        real(dp) :: norm

        ierr = ERROR_SUCCESS

        if (n_deg <= 1) then
            ierr = ERROR_SUCCESS ! Nothing to do for a single vector
            return
        end if

        if (L <= 0) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (size(indexes) < n_deg) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        if (size(eigvecs, 1) /= L .or. size(eigvecs, 2) < L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        do k = 1, n_deg
            if (indexes(k) < 1 .or. indexes(k) > L) then
                ierr = ERROR_INVALID_INPUT
                return
            end if
        end do

        ! Modified Gram-Schmidt orthonormalization
        do k = 1, n_deg
            idx_k = indexes(k)

            do j = 1, k - 1
                idx_j = indexes(j)
                proj = DOT_PRODUCT(eigvecs(:,idx_j), eigvecs(:,idx_k))
                eigvecs(:,idx_k) = eigvecs(:,idx_k) - proj * eigvecs(:,idx_j)
            end do

            norm = SQRT(real(DOT_PRODUCT(eigvecs(:,idx_k), eigvecs(:,idx_k)), kind=dp))
            if (norm < 1.0e-14_dp) then
                ierr = ERROR_LINEAR_DEPENDENCE
                return
            end if

            eigvecs(:,idx_k) = eigvecs(:,idx_k) / norm
        end do
    end subroutine orthonormalize_degenerate_subspace_complex

    !> @brief Count number of eigenvalues degenerate with given index
    !!
    !! Searches forward and backward from given index to count how many
    !! eigenvalues are within tolerance.
    !!
    !! @param[in] eigvals Eigenvalues (sorted)
    !! @param[in] L Number of eigenvalues
    !! @param[in] index Index of reference eigenvalue
    !! @return Number of degenerate eigenvalues (including reference)
    function compute_degeneracy_count(eigvals, L, index) result(deg_count)
        real(dp), intent(in) :: eigvals(:)
        integer, intent(in) :: L, index
        real(dp) :: eval
        integer :: deg_count, i

        if (index < 1 .or. index > L .or. L <= 0) then
            deg_count = 0
            return
        end if

        eval = eigvals(index)
        deg_count = 1

        do i = index - 1, 1, -1
            if (abs(eigvals(i) - eval) < DEG_TOL) then
                deg_count = deg_count + 1
            else
                exit
            end if
        end do

        do i = index + 1, L
            if (abs(eigvals(i) - eval) < DEG_TOL) then
                deg_count = deg_count + 1
            else
                exit
            end if
        end do
    end function compute_degeneracy_count

    !> @brief Verify eigenvectors are orthonormal
    !!
    !! Computes overlap matrix S = V^T * V using BLAS and checks
    !! deviation from identity matrix.
    !!
    !! @param[in] eigvecs Eigenvector matrix (L × L)
    !! @param[in] L Dimension
    !! @param[out] is_orthonormal True if ||S - I||_∞ < tol
    !! @param[out] max_deviation Maximum deviation from orthonormality
    !! @param[out] ierr Error code
    subroutine verify_orthonormality(eigvecs, L, is_orthonormal, &
                                  max_deviation, ierr)
        real(dp), intent(in) :: eigvecs(:,:)
        integer, intent(in) :: L
        logical, intent(out) :: is_orthonormal
        real(dp), intent(out) :: max_deviation
        integer, intent(out) :: ierr
        
        real(dp), allocatable :: S(:,:)
        real(dp) :: dev
        integer :: i, j
        
        ierr = ERROR_SUCCESS

        if (L <= 0) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (size(eigvecs, 1) /= L .or. size(eigvecs, 2) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        allocate(S(L, L))
        call DGEMM('T',        & ! Transpose first matrix
                   'N',        & ! Don't transpose second
                   L, L, L,    & ! Dimensions: m, n, k
                   1.0_dp,     & ! alpha
                   eigvecs, L, & ! A matrix
                   eigvecs, L, & ! B matrix
                   0.0_dp,     & ! beta
                   S, L)         ! C = result                        
        
        max_deviation = 0.0_dp
        do i = 1, L
            do j = 1, L
                if (i == j) then
                    dev = abs(S(i,j) - 1.0_dp)
                else
                    dev = abs(S(i,j))
                end if
                max_deviation = max(max_deviation, dev)
            end do
        end do
        
        is_orthonormal = (max_deviation < 1.0e-8_dp)
        deallocate(S)
    end subroutine verify_orthonormality
end module degeneracy_handler