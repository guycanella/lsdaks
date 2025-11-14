!> Unit tests for degeneracy_handler module
program test_degeneracy_handler
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp
    implicit none

    real(dp), parameter :: TOL = 1.0e-9_dp

    call execute_serial_cmd_app(get_degeneracy_tests())

contains

    function get_degeneracy_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests

        tests = test_list([ &
            test("find_subspaces_no_degeneracy", test_find_subspaces_no_degeneracy), &
            test("find_subspaces_single_pair", test_find_subspaces_single_pair), &
            test("find_subspaces_triple_degeneracy", test_find_subspaces_triple_degeneracy), &
            test("find_subspaces_multiple_groups", test_find_subspaces_multiple_groups), &
            test("find_subspaces_all_degenerate", test_find_subspaces_all_degenerate), &
            test("degeneracy_count_single", test_degeneracy_count_single), &
            test("degeneracy_count_pair", test_degeneracy_count_pair), &
            test("degeneracy_count_triple", test_degeneracy_count_triple), &
            test("orthonormalize_real_pair", test_orthonormalize_real_pair), &
            test("orthonormalize_complex_pair", test_orthonormalize_complex_pair), &
            test("verify_orthonormality_identity", test_verify_orthonormality_identity), &
            test("verify_orthonormality_non_orthogonal", test_verify_orthonormality_non_orthogonal) &
        ])
    end function get_degeneracy_tests

    !> Test finding degenerate subspaces with no degeneracies
    !!
    !! Physics: For a generic Hamiltonian without symmetries, eigenvalues
    !! are typically non-degenerate. Each energy level is unique.
    subroutine test_find_subspaces_no_degeneracy()
        use fortuno_serial, only: check => serial_check
        use degeneracy_handler, only: find_degenerate_subspaces
        integer, parameter :: L = 5
        real(dp) :: eigvals(L)
        integer, allocatable :: subspaces(:,:), subspace_sizes(:)
        integer :: n_subspaces, ierr

        ! All distinct eigenvalues
        eigvals = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]

        call find_degenerate_subspaces(eigvals, L, subspaces, subspace_sizes, &
                                       n_subspaces, ierr)

        call check(ierr == 0, "Should succeed with no degeneracies")
        call check(n_subspaces == 0, "Should find no degenerate subspaces")
        call check(size(subspaces, 1) == 0, "Subspaces array should be empty")
    end subroutine test_find_subspaces_no_degeneracy

    !> Test finding single pair of degenerate eigenvalues
    !!
    !! Physics: Two-fold degeneracy occurs when a system has a discrete
    !! symmetry (e.g., spin up/down, parity even/odd). Common in quantum
    !! systems.
    subroutine test_find_subspaces_single_pair()
        use fortuno_serial, only: check => serial_check
        use degeneracy_handler, only: find_degenerate_subspaces
        integer, parameter :: L = 5
        real(dp) :: eigvals(L)
        integer, allocatable :: subspaces(:,:), subspace_sizes(:)
        integer :: n_subspaces, ierr

        ! Two degenerate eigenvalues at positions 2 and 3
        eigvals = [1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]

        call find_degenerate_subspaces(eigvals, L, subspaces, subspace_sizes, &
                                       n_subspaces, ierr)

        call check(ierr == 0, "Should succeed")
        call check(n_subspaces == 1, "Should find one degenerate subspace")
        call check(subspace_sizes(1) == 2, "Subspace should have size 2")
        call check(subspaces(1,1) == 2, "First index should be 2")
        call check(subspaces(1,2) == 3, "Second index should be 3")
    end subroutine test_find_subspaces_single_pair

    !> Test finding triple degeneracy
    !!
    !! Physics: Three-fold degeneracy can occur in systems with higher
    !! symmetries, like angular momentum l=1 (px, py, pz orbitals) or
    !! cubic lattices.
    subroutine test_find_subspaces_triple_degeneracy()
        use fortuno_serial, only: check => serial_check
        use degeneracy_handler, only: find_degenerate_subspaces
        integer, parameter :: L = 6
        real(dp) :: eigvals(L)
        integer, allocatable :: subspaces(:,:), subspace_sizes(:)
        integer :: n_subspaces, ierr

        ! Three degenerate eigenvalues at positions 3, 4, 5
        eigvals = [1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 4.0_dp]

        call find_degenerate_subspaces(eigvals, L, subspaces, subspace_sizes, &
                                       n_subspaces, ierr)

        call check(ierr == 0, "Should succeed")
        call check(n_subspaces == 1, "Should find one degenerate subspace")
        call check(subspace_sizes(1) == 3, "Subspace should have size 3")
        call check(subspaces(1,1) == 3, "First index should be 3")
        call check(subspaces(1,2) == 4, "Second index should be 4")
        call check(subspaces(1,3) == 5, "Third index should be 5")
    end subroutine test_find_subspaces_triple_degeneracy

    !> Test finding multiple degenerate groups
    !!
    !! Physics: Complex systems can have multiple independent degenerate
    !! subspaces at different energies (e.g., different angular momentum
    !! manifolds in atoms).
    subroutine test_find_subspaces_multiple_groups()
        use fortuno_serial, only: check => serial_check
        use degeneracy_handler, only: find_degenerate_subspaces
        integer, parameter :: L = 8
        real(dp) :: eigvals(L)
        integer, allocatable :: subspaces(:,:), subspace_sizes(:)
        integer :: n_subspaces, ierr

        ! Two pairs of degeneracies: (1,2) and (5,6)
        eigvals = [1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]

        call find_degenerate_subspaces(eigvals, L, subspaces, subspace_sizes, &
                                       n_subspaces, ierr)

        call check(ierr == 0, "Should succeed")
        call check(n_subspaces == 2, "Should find two degenerate subspaces")
        call check(subspace_sizes(1) == 2, "First subspace size 2")
        call check(subspace_sizes(2) == 2, "Second subspace size 2")
        call check(subspaces(1,1) == 1 .and. subspaces(1,2) == 2, "First pair at (1,2)")
        call check(subspaces(2,1) == 5 .and. subspaces(2,2) == 6, "Second pair at (5,6)")
    end subroutine test_find_subspaces_multiple_groups

    !> Test all eigenvalues degenerate
    !!
    !! Physics: Maximum degeneracy occurs for the identity operator,
    !! where all states have the same energy. Rare but mathematically
    !! valid case.
    subroutine test_find_subspaces_all_degenerate()
        use fortuno_serial, only: check => serial_check
        use degeneracy_handler, only: find_degenerate_subspaces
        integer, parameter :: L = 4
        real(dp) :: eigvals(L)
        integer, allocatable :: subspaces(:,:), subspace_sizes(:)
        integer :: n_subspaces, ierr

        eigvals = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]

        call find_degenerate_subspaces(eigvals, L, subspaces, subspace_sizes, &
                                       n_subspaces, ierr)

        call check(ierr == 0, "Should succeed")
        call check(n_subspaces == 1, "Should find one subspace")
        call check(subspace_sizes(1) == 4, "All 4 eigenvalues degenerate")
    end subroutine test_find_subspaces_all_degenerate

    !> Test degeneracy count for non-degenerate eigenvalue
    !!
    !! Physics: For a unique energy level, degeneracy count is 1.
    subroutine test_degeneracy_count_single()
        use fortuno_serial, only: check => serial_check
        use degeneracy_handler, only: compute_degeneracy_count
        integer, parameter :: L = 5
        real(dp) :: eigvals(L)
        integer :: deg_count

        eigvals = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]

        deg_count = compute_degeneracy_count(eigvals, L, 3)
        call check(deg_count == 1, "Non-degenerate eigenvalue should have count 1")
    end subroutine test_degeneracy_count_single

    !> Test degeneracy count for pair
    !!
    !! Physics: Two-fold degenerate level has count 2.
    subroutine test_degeneracy_count_pair()
        use fortuno_serial, only: check => serial_check
        use degeneracy_handler, only: compute_degeneracy_count
        integer, parameter :: L = 5
        real(dp) :: eigvals(L)
        integer :: deg_count

        eigvals = [1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]

        deg_count = compute_degeneracy_count(eigvals, L, 2)
        call check(deg_count == 2, "Pair degeneracy from index 2")

        deg_count = compute_degeneracy_count(eigvals, L, 3)
        call check(deg_count == 2, "Pair degeneracy from index 3")
    end subroutine test_degeneracy_count_pair

    !> Test degeneracy count for triple
    !!
    !! Physics: Three-fold degenerate level has count 3 from any index.
    subroutine test_degeneracy_count_triple()
        use fortuno_serial, only: check => serial_check
        use degeneracy_handler, only: compute_degeneracy_count
        integer, parameter :: L = 6
        real(dp) :: eigvals(L)
        integer :: deg_count

        eigvals = [1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 4.0_dp]

        deg_count = compute_degeneracy_count(eigvals, L, 3)
        call check(deg_count == 3, "Triple degeneracy from index 3")

        deg_count = compute_degeneracy_count(eigvals, L, 4)
        call check(deg_count == 3, "Triple degeneracy from index 4")

        deg_count = compute_degeneracy_count(eigvals, L, 5)
        call check(deg_count == 3, "Triple degeneracy from index 5")
    end subroutine test_degeneracy_count_triple

    !> Test real eigenvector re-orthonormalization for pair
    !!
    !! Physics: When two eigenstates are degenerate, LAPACK may return
    !! an arbitrary basis. QR decomposition provides a stable way to
    !! obtain an orthonormal basis for the degenerate subspace.
    subroutine test_orthonormalize_real_pair()
        use fortuno_serial, only: check => serial_check
        use degeneracy_handler, only: orthonormalize_degenerate_subspace
        integer, parameter :: L = 4
        real(dp) :: eigvecs(L,L)
        integer :: indices(2)
        real(dp) :: overlap
        integer :: ierr

        ! Create two non-orthogonal vectors at indices 2 and 3
        eigvecs = 0.0_dp
        eigvecs(1,2) = 1.0_dp
        eigvecs(2,2) = 1.0_dp
        eigvecs(1,3) = 1.0_dp
        eigvecs(2,3) = 0.5_dp

        indices = [2, 3]

        call orthonormalize_degenerate_subspace(eigvecs, L, indices, 2, ierr)

        call check(ierr == 0, "Orthonormalization should succeed")

        ! Check normalization
        call check(abs(sqrt(dot_product(eigvecs(:,2), eigvecs(:,2))) - 1.0_dp) < TOL, &
                   "First vector should be normalized")
        call check(abs(sqrt(dot_product(eigvecs(:,3), eigvecs(:,3))) - 1.0_dp) < TOL, &
                   "Second vector should be normalized")

        ! Check orthogonality
        overlap = abs(dot_product(eigvecs(:,2), eigvecs(:,3)))
        call check(overlap < TOL, "Vectors should be orthogonal")
    end subroutine test_orthonormalize_real_pair

    !> Test complex eigenvector re-orthonormalization
    !!
    !! Physics: For twisted boundary conditions or spin-orbit coupling,
    !! eigenvectors are complex. Modified Gram-Schmidt provides orthonormal
    !! basis in the degenerate subspace.
    subroutine test_orthonormalize_complex_pair()
        use fortuno_serial, only: check => serial_check
        use degeneracy_handler, only: orthonormalize_degenerate_subspace_complex
        integer, parameter :: L = 4
        complex(dp) :: eigvecs(L,L)
        integer :: indices(2)
        complex(dp) :: overlap
        real(dp) :: norm
        integer :: ierr

        ! Create two non-orthogonal complex vectors
        eigvecs = cmplx(0.0_dp, 0.0_dp, kind=dp)
        eigvecs(1,2) = cmplx(1.0_dp, 0.0_dp, kind=dp)
        eigvecs(2,2) = cmplx(0.0_dp, 1.0_dp, kind=dp)
        eigvecs(1,3) = cmplx(1.0_dp, 1.0_dp, kind=dp)
        eigvecs(2,3) = cmplx(0.5_dp, 0.0_dp, kind=dp)

        indices = [2, 3]

        call orthonormalize_degenerate_subspace_complex(eigvecs, L, indices, 2, ierr)

        call check(ierr == 0, "Complex orthonormalization should succeed")

        ! Check normalization
        norm = sqrt(real(dot_product(eigvecs(:,2), eigvecs(:,2)), kind=dp))
        call check(abs(norm - 1.0_dp) < TOL, "First vector normalized")

        norm = sqrt(real(dot_product(eigvecs(:,3), eigvecs(:,3)), kind=dp))
        call check(abs(norm - 1.0_dp) < TOL, "Second vector normalized")

        ! Check orthogonality
        overlap = dot_product(eigvecs(:,2), eigvecs(:,3))
        call check(abs(overlap) < TOL, "Vectors should be orthogonal")
    end subroutine test_orthonormalize_complex_pair

    !> Test verification of orthonormality for identity matrix
    !!
    !! Physics: The standard basis (identity matrix) is perfectly
    !! orthonormal. This should be detected with zero deviation.
    subroutine test_verify_orthonormality_identity()
        use fortuno_serial, only: check => serial_check
        use degeneracy_handler, only: verify_orthonormality
        integer, parameter :: L = 5
        real(dp) :: eigvecs(L,L)
        logical :: is_orthonormal
        real(dp) :: max_deviation
        integer :: ierr, i

        ! Identity matrix (standard basis)
        eigvecs = 0.0_dp
        do i = 1, L
            eigvecs(i,i) = 1.0_dp
        end do

        call verify_orthonormality(eigvecs, L, is_orthonormal, max_deviation, ierr)

        call check(ierr == 0, "Verification should succeed")
        call check(is_orthonormal, "Identity should be orthonormal")
        call check(max_deviation < TOL, "Deviation should be near zero")
    end subroutine test_verify_orthonormality_identity

    !> Test verification detects non-orthogonal vectors
    !!
    !! Physics: If eigenvectors are not properly orthogonalized, physical
    !! calculations (overlaps, expectation values) will be incorrect.
    !! This test ensures we can detect such problems.
    subroutine test_verify_orthonormality_non_orthogonal()
        use fortuno_serial, only: check => serial_check
        use degeneracy_handler, only: verify_orthonormality
        integer, parameter :: L = 3
        real(dp) :: eigvecs(L,L)
        logical :: is_orthonormal
        real(dp) :: max_deviation
        integer :: ierr

        ! Create non-orthogonal but normalized vectors
        eigvecs = 0.0_dp
        eigvecs(1,1) = 1.0_dp
        eigvecs(1,2) = 0.6_dp
        eigvecs(2,2) = 0.8_dp
        eigvecs(3,3) = 1.0_dp

        call verify_orthonormality(eigvecs, L, is_orthonormal, max_deviation, ierr)

        call check(ierr == 0, "Verification should succeed")
        call check(.not. is_orthonormal, "Should detect non-orthogonality")
        call check(max_deviation > TOL, "Deviation should be significant")
    end subroutine test_verify_orthonormality_non_orthogonal

end program test_degeneracy_handler
