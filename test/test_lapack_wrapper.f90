!> Unit tests for lapack_wrapper module
program test_lapack_wrapper
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp, PI
    implicit none

    real(dp), parameter :: TOL = 1.0e-9_dp

    call execute_serial_cmd_app(get_lapack_tests())

contains

    function get_lapack_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests

        tests = test_list([ &
            test("validate_inputs_valid", test_validate_inputs_valid), &
            test("validate_inputs_invalid_L", test_validate_inputs_invalid_L), &
            test("validate_inputs_size_mismatch_H", test_validate_inputs_size_mismatch_H), &
            test("validate_inputs_size_mismatch_eigvals", test_validate_inputs_size_mismatch_eigvals), &
            test("validate_inputs_size_mismatch_eigvecs", test_validate_inputs_size_mismatch_eigvecs), &
            test("validate_inputs_eigvecs_not_present", test_validate_inputs_eigvecs_not_present), &
            test("diag_real_identity_matrix", test_diag_real_identity_matrix), &
            test("diag_real_diagonal_matrix", test_diag_real_diagonal_matrix), &
            test("diag_real_symmetric_2x2", test_diag_real_symmetric_2x2), &
            test("diag_real_tridiagonal", test_diag_real_tridiagonal), &
            test("diag_real_eigenvalue_order", test_diag_real_eigenvalue_order), &
            test("diag_real_eigenvector_normalization", test_diag_real_eigenvector_normalization), &
            test("diag_real_eigenvector_orthogonality", test_diag_real_eigenvector_orthogonality), &
            test("diag_real_values_only", test_diag_real_values_only), &
            test("diag_complex_identity", test_diag_complex_identity), &
            test("diag_complex_hermitian", test_diag_complex_hermitian), &
            test("diag_complex_eigenvalues_real", test_diag_complex_eigenvalues_real), &
            test("diag_complex_values_only", test_diag_complex_values_only) &
        ])
    end function get_lapack_tests

    !> Test validation with correct inputs
    !!
    !! Physics: Before diagonalization, verify all input arrays have compatible
    !! dimensions to prevent memory access errors and ensure physically meaningful
    !! results.
    subroutine test_validate_inputs_valid()
        use fortuno_serial, only: check => serial_check
        use lapack_wrapper, only: validate_diagonalization_inputs
        integer, parameter :: L = 10
        real(dp) :: H(L,L), eigvals(L), eigvecs(L,L)
        integer :: ierr

        H = 0.0_dp
        call validate_diagonalization_inputs(L, H, eigvals, eigvecs, .true., ierr)
        call check(ierr == 0, "Valid inputs with eigenvectors should pass")

        call validate_diagonalization_inputs(L, H, eigvals, compute_vectors=.false., ierr=ierr)
        call check(ierr == 0, "Valid inputs without eigenvectors should pass")
    end subroutine test_validate_inputs_valid

    !> Test validation failure for invalid L
    !!
    !! Physics: A Hamiltonian must represent at least one lattice site (L > 0).
    !! L ≤ 0 is unphysical and indicates a programming error.
    subroutine test_validate_inputs_invalid_L()
        use fortuno_serial, only: check => serial_check
        use lapack_wrapper, only: validate_diagonalization_inputs
        use lsda_errors, only: ERROR_INVALID_INPUT
        real(dp) :: H(1,1), eigvals(1), eigvecs(1,1)
        integer :: ierr

        call validate_diagonalization_inputs(0, H, eigvals, eigvecs, .true., ierr)
        call check(ierr == ERROR_INVALID_INPUT, "L=0 should fail")

        call validate_diagonalization_inputs(-5, H, eigvals, eigvecs, .true., ierr)
        call check(ierr == ERROR_INVALID_INPUT, "L<0 should fail")
    end subroutine test_validate_inputs_invalid_L

    !> Test validation failure for H size mismatch
    !!
    !! Physics: The Hamiltonian matrix must be square (L×L) to represent
    !! the tight-binding model on L sites properly.
    subroutine test_validate_inputs_size_mismatch_H()
        use fortuno_serial, only: check => serial_check
        use lapack_wrapper, only: validate_diagonalization_inputs
        use lsda_errors, only: ERROR_SIZE_MISMATCH
        integer, parameter :: L = 5
        real(dp) :: H(L+1,L), eigvals(L), eigvecs(L,L)
        integer :: ierr

        call validate_diagonalization_inputs(L, H, eigvals, eigvecs, .true., ierr)
        call check(ierr == ERROR_SIZE_MISMATCH, "H size mismatch should fail")
    end subroutine test_validate_inputs_size_mismatch_H

    !> Test validation failure for eigenvalue array size mismatch
    !!
    !! Physics: Must have exactly L eigenvalues for an L×L matrix.
    subroutine test_validate_inputs_size_mismatch_eigvals()
        use fortuno_serial, only: check => serial_check
        use lapack_wrapper, only: validate_diagonalization_inputs
        use lsda_errors, only: ERROR_SIZE_MISMATCH
        integer, parameter :: L = 5
        real(dp) :: H(L,L), eigvals(L-1), eigvecs(L,L)
        integer :: ierr

        call validate_diagonalization_inputs(L, H, eigvals, eigvecs, .true., ierr)
        call check(ierr == ERROR_SIZE_MISMATCH, "eigvals size mismatch should fail")
    end subroutine test_validate_inputs_size_mismatch_eigvals

    !> Test validation failure for eigenvector array size mismatch
    !!
    !! Physics: Eigenvector matrix must be L×L to store all L eigenvectors
    !! of dimension L.
    subroutine test_validate_inputs_size_mismatch_eigvecs()
        use fortuno_serial, only: check => serial_check
        use lapack_wrapper, only: validate_diagonalization_inputs
        use lsda_errors, only: ERROR_SIZE_MISMATCH
        integer, parameter :: L = 5
        real(dp) :: H(L,L), eigvals(L), eigvecs(L,L-1)
        integer :: ierr

        call validate_diagonalization_inputs(L, H, eigvals, eigvecs, .true., ierr)
        call check(ierr == ERROR_SIZE_MISMATCH, "eigvecs size mismatch should fail")
    end subroutine test_validate_inputs_size_mismatch_eigvecs

    !> Test validation failure when eigenvectors requested but not provided
    !!
    !! Physics: If eigenvectors are needed (compute_vectors=.true.), the
    !! eigvecs array must be present to store the results.
    subroutine test_validate_inputs_eigvecs_not_present()
        use fortuno_serial, only: check => serial_check
        use lapack_wrapper, only: validate_diagonalization_inputs
        use lsda_errors, only: ERROR_INVALID_INPUT
        integer, parameter :: L = 5
        real(dp) :: H(L,L), eigvals(L)
        integer :: ierr

        call validate_diagonalization_inputs(L, H, eigvals, compute_vectors=.true., ierr=ierr)
        call check(ierr == ERROR_INVALID_INPUT, "Missing eigvecs should fail when compute_vectors=.true.")
    end subroutine test_validate_inputs_eigvecs_not_present

    !> Test diagonalization of identity matrix
    !!
    !! Physics: The identity matrix I has eigenvalues all equal to 1,
    !! with any orthonormal basis as eigenvectors. This is the trivial
    !! case where all states are degenerate.
    subroutine test_diag_real_identity_matrix()
        use fortuno_serial, only: check => serial_check
        use lapack_wrapper, only: diagonalize_symmetric_real
        integer, parameter :: L = 4
        real(dp) :: H(L,L), eigvals(L), eigvecs(L,L)
        integer :: ierr, i

        H = 0.0_dp
        do i = 1, L
            H(i,i) = 1.0_dp
        end do

        call diagonalize_symmetric_real(H, L, eigvals, eigvecs, ierr)

        call check(ierr == 0, "Diagonalization should succeed")
        do i = 1, L
            call check(abs(eigvals(i) - 1.0_dp) < TOL, "All eigenvalues should be 1")
        end do
    end subroutine test_diag_real_identity_matrix

    !> Test diagonalization of diagonal matrix
    !!
    !! Physics: A diagonal matrix already has eigenvalues on the diagonal
    !! and eigenvectors as standard basis vectors. Tests that LAPACK
    !! correctly identifies this trivial case.
    subroutine test_diag_real_diagonal_matrix()
        use fortuno_serial, only: check => serial_check
        use lapack_wrapper, only: diagonalize_symmetric_real
        integer, parameter :: L = 5
        real(dp) :: H(L,L), eigvals(L), eigvecs(L,L)
        real(dp) :: diag_values(L)
        integer :: ierr, i

        diag_values = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]

        H = 0.0_dp
        do i = 1, L
            H(i,i) = diag_values(i)
        end do

        call diagonalize_symmetric_real(H, L, eigvals, eigvecs, ierr)

        call check(ierr == 0, "Diagonalization should succeed")
        ! Eigenvalues should be sorted in ascending order
        do i = 1, L
            call check(abs(eigvals(i) - real(i, dp)) < TOL, "Eigenvalues should match diagonal")
        end do
    end subroutine test_diag_real_diagonal_matrix

    !> Test diagonalization of 2×2 symmetric matrix with known eigenvalues
    !!
    !! Physics: For a 2×2 symmetric matrix, eigenvalues can be computed
    !! analytically. This validates numerical diagonalization against
    !! exact solution.
    subroutine test_diag_real_symmetric_2x2()
        use fortuno_serial, only: check => serial_check
        use lapack_wrapper, only: diagonalize_symmetric_real
        integer, parameter :: L = 2
        real(dp) :: H(L,L), eigvals(L), eigvecs(L,L)
        real(dp) :: lambda1, lambda2
        integer :: ierr

        ! H = [1  2]
        !     [2  1]
        ! Eigenvalues: 3, -1
        H(1,1) = 1.0_dp
        H(1,2) = 2.0_dp
        H(2,1) = 2.0_dp
        H(2,2) = 1.0_dp

        call diagonalize_symmetric_real(H, L, eigvals, eigvecs, ierr)

        call check(ierr == 0, "Diagonalization should succeed")
        call check(abs(eigvals(1) + 1.0_dp) < TOL, "First eigenvalue should be -1")
        call check(abs(eigvals(2) - 3.0_dp) < TOL, "Second eigenvalue should be 3")
    end subroutine test_diag_real_symmetric_2x2

    !> Test diagonalization of tridiagonal tight-binding Hamiltonian
    !!
    !! Physics: The tight-binding Hamiltonian with open BC is tridiagonal:
    !! H(i,i) = 0, H(i,i±1) = -1. Eigenvalues are E_n = -2cos(nπ/(L+1)).
    !! This validates against known analytical solution.
    subroutine test_diag_real_tridiagonal()
        use fortuno_serial, only: check => serial_check
        use lapack_wrapper, only: diagonalize_symmetric_real
        integer, parameter :: L = 10
        real(dp) :: H(L,L), eigvals(L), eigvecs(L,L)
        real(dp) :: E_exact
        integer :: ierr, i

        ! Build tight-binding Hamiltonian
        H = 0.0_dp
        do i = 1, L-1
            H(i,i+1) = -1.0_dp
            H(i+1,i) = -1.0_dp
        end do

        call diagonalize_symmetric_real(H, L, eigvals, eigvecs, ierr)

        call check(ierr == 0, "Diagonalization should succeed")

        ! Check ground state eigenvalue (n=1)
        E_exact = -2.0_dp * cos(PI / real(L+1, dp))
        call check(abs(eigvals(1) - E_exact) < TOL, &
                   "Ground state eigenvalue should match analytical formula")

        ! All eigenvalues must be in band [-2, 2]
        call check(all(eigvals >= -2.0_dp - TOL), "All eigenvalues >= -2")
        call check(all(eigvals <= 2.0_dp + TOL), "All eigenvalues <= 2")
    end subroutine test_diag_real_tridiagonal

    !> Test that eigenvalues are returned in ascending order
    !!
    !! Physics: Standard convention is to order eigenvalues from lowest
    !! to highest energy. Ground state is eigenvalue #1.
    subroutine test_diag_real_eigenvalue_order()
        use fortuno_serial, only: check => serial_check
        use lapack_wrapper, only: diagonalize_symmetric_real
        integer, parameter :: L = 6
        real(dp) :: H(L,L), eigvals(L), eigvecs(L,L)
        integer :: ierr, i

        ! Random-looking diagonal matrix
        H = 0.0_dp
        H(1,1) = 5.0_dp
        H(2,2) = 1.0_dp
        H(3,3) = 3.0_dp
        H(4,4) = 2.0_dp
        H(5,5) = 6.0_dp
        H(6,6) = 4.0_dp

        call diagonalize_symmetric_real(H, L, eigvals, eigvecs, ierr)

        call check(ierr == 0, "Diagonalization should succeed")

        ! Check ascending order
        do i = 1, L-1
            call check(eigvals(i) <= eigvals(i+1) + TOL, "Eigenvalues should be sorted ascending")
        end do
    end subroutine test_diag_real_eigenvalue_order

    !> Test eigenvector normalization
    !!
    !! Physics: Eigenvectors must be normalized (||v|| = 1) to form
    !! a proper orthonormal basis for the Hilbert space.
    subroutine test_diag_real_eigenvector_normalization()
        use fortuno_serial, only: check => serial_check
        use lapack_wrapper, only: diagonalize_symmetric_real
        integer, parameter :: L = 5
        real(dp) :: H(L,L), eigvals(L), eigvecs(L,L)
        real(dp) :: norm
        integer :: ierr, i

        H = 0.0_dp
        do i = 1, L-1
            H(i,i+1) = -1.0_dp
            H(i+1,i) = -1.0_dp
        end do

        call diagonalize_symmetric_real(H, L, eigvals, eigvecs, ierr)

        call check(ierr == 0, "Diagonalization should succeed")

        ! Check normalization: v^T * v = 1
        do i = 1, L
            norm = sqrt(dot_product(eigvecs(:,i), eigvecs(:,i)))
            call check(abs(norm - 1.0_dp) < TOL, "Eigenvectors should be normalized")
        end do
    end subroutine test_diag_real_eigenvector_normalization

    !> Test eigenvector orthogonality
    !!
    !! Physics: Eigenvectors corresponding to different eigenvalues
    !! of a Hermitian operator are orthogonal. This is fundamental
    !! to quantum mechanics.
    subroutine test_diag_real_eigenvector_orthogonality()
        use fortuno_serial, only: check => serial_check
        use lapack_wrapper, only: diagonalize_symmetric_real
        integer, parameter :: L = 5
        real(dp) :: H(L,L), eigvals(L), eigvecs(L,L)
        real(dp) :: overlap
        integer :: ierr, i, j

        H = 0.0_dp
        do i = 1, L-1
            H(i,i+1) = -1.0_dp
            H(i+1,i) = -1.0_dp
        end do

        call diagonalize_symmetric_real(H, L, eigvals, eigvecs, ierr)

        call check(ierr == 0, "Diagonalization should succeed")

        ! Check orthogonality: v_i^T * v_j = 0 for i ≠ j
        do i = 1, L
            do j = i+1, L
                overlap = abs(dot_product(eigvecs(:,i), eigvecs(:,j)))
                call check(overlap < TOL, "Eigenvectors should be orthogonal")
            end do
        end do
    end subroutine test_diag_real_eigenvector_orthogonality

    !> Test eigenvalues-only version is faster and correct
    !!
    !! Physics: When only energy levels are needed (not wavefunctions),
    !! we can skip eigenvector computation for ~2x speedup.
    subroutine test_diag_real_values_only()
        use fortuno_serial, only: check => serial_check
        use lapack_wrapper, only: diagonalize_symmetric_real_values_only
        integer, parameter :: L = 8
        real(dp) :: H(L,L), eigvals(L)
        real(dp) :: E_exact
        integer :: ierr, i

        H = 0.0_dp
        do i = 1, L-1
            H(i,i+1) = -1.0_dp
            H(i+1,i) = -1.0_dp
        end do

        call diagonalize_symmetric_real_values_only(H, L, eigvals, ierr)

        call check(ierr == 0, "Eigenvalues-only diagonalization should succeed")

        ! Check ground state
        E_exact = -2.0_dp * cos(PI / real(L+1, dp))
        call check(abs(eigvals(1) - E_exact) < TOL, &
                   "Ground state energy should match")
    end subroutine test_diag_real_values_only

    !> Test diagonalization of complex Hermitian identity matrix
    !!
    !! Physics: Even for complex matrices, a Hermitian matrix (H† = H)
    !! has real eigenvalues. The identity is the simplest Hermitian matrix.
    subroutine test_diag_complex_identity()
        use fortuno_serial, only: check => serial_check
        use lapack_wrapper, only: diagonalize_hermitian_complex
        integer, parameter :: L = 4
        complex(dp) :: H(L,L), eigvecs(L,L)
        real(dp) :: eigvals(L)
        integer :: ierr, i

        H = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do i = 1, L
            H(i,i) = cmplx(1.0_dp, 0.0_dp, kind=dp)
        end do

        call diagonalize_hermitian_complex(H, L, eigvals, eigvecs, ierr)

        call check(ierr == 0, "Complex diagonalization should succeed")
        do i = 1, L
            call check(abs(eigvals(i) - 1.0_dp) < TOL, "Eigenvalues should be 1")
        end do
    end subroutine test_diag_complex_identity

    !> Test diagonalization of complex Hermitian matrix
    !!
    !! Physics: For twisted boundary conditions with Aharonov-Bohm phase,
    !! the Hamiltonian becomes complex. However, H† = H ensures real
    !! eigenvalues (energies).
    subroutine test_diag_complex_hermitian()
        use fortuno_serial, only: check => serial_check
        use lapack_wrapper, only: diagonalize_hermitian_complex
        integer, parameter :: L = 2
        complex(dp) :: H(L,L), eigvecs(L,L)
        real(dp) :: eigvals(L)
        integer :: ierr

        ! H = [1    i ]
        !     [-i   1 ]
        ! Hermitian, eigenvalues: 0, 2
        H(1,1) = cmplx(1.0_dp, 0.0_dp, kind=dp)
        H(1,2) = cmplx(0.0_dp, 1.0_dp, kind=dp)
        H(2,1) = cmplx(0.0_dp, -1.0_dp, kind=dp)
        H(2,2) = cmplx(1.0_dp, 0.0_dp, kind=dp)

        call diagonalize_hermitian_complex(H, L, eigvals, eigvecs, ierr)

        call check(ierr == 0, "Complex diagonalization should succeed")
        call check(abs(eigvals(1) - 0.0_dp) < TOL, "First eigenvalue should be 0")
        call check(abs(eigvals(2) - 2.0_dp) < TOL, "Second eigenvalue should be 2")
    end subroutine test_diag_complex_hermitian

    !> Test that complex Hermitian eigenvalues are real
    !!
    !! Physics: A fundamental theorem of quantum mechanics states that
    !! observables (Hermitian operators) have real eigenvalues (measurable
    !! energies). This must hold even when H is complex.
    subroutine test_diag_complex_eigenvalues_real()
        use fortuno_serial, only: check => serial_check
        use lapack_wrapper, only: diagonalize_hermitian_complex
        integer, parameter :: L = 5
        complex(dp) :: H(L,L), eigvecs(L,L)
        real(dp) :: eigvals(L)
        integer :: ierr, i

        ! Build tight-binding with phase (twisted BC)
        H = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do i = 1, L-1
            H(i,i+1) = cmplx(-1.0_dp, 0.0_dp, kind=dp)
            H(i+1,i) = cmplx(-1.0_dp, 0.0_dp, kind=dp)
        end do
        ! Add twist: H(1,L) = -exp(iθ), H(L,1) = -exp(-iθ)
        H(1,L) = cmplx(-0.707_dp, -0.707_dp, kind=dp)  ! θ = π/4
        H(L,1) = cmplx(-0.707_dp, 0.707_dp, kind=dp)

        call diagonalize_hermitian_complex(H, L, eigvals, eigvecs, ierr)

        call check(ierr == 0, "Twisted BC diagonalization should succeed")
        ! All eigenvalues should be in band [-2, 2]
        call check(all(eigvals >= -2.0_dp - TOL), "Eigenvalues >= -2")
        call check(all(eigvals <= 2.0_dp + TOL), "Eigenvalues <= 2")
    end subroutine test_diag_complex_eigenvalues_real

    !> Test complex eigenvalues-only version
    !!
    !! Physics: For complex Hamiltonians, eigenvalues are still real
    !! (if Hermitian). Values-only version should match full diagonalization.
    subroutine test_diag_complex_values_only()
        use fortuno_serial, only: check => serial_check
        use lapack_wrapper, only: diagonalize_hermitian_complex_values_only
        integer, parameter :: L = 6
        complex(dp) :: H(L,L)
        real(dp) :: eigvals(L)
        integer :: ierr, i

        H = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do i = 1, L-1
            H(i,i+1) = cmplx(-1.0_dp, 0.0_dp, kind=dp)
            H(i+1,i) = cmplx(-1.0_dp, 0.0_dp, kind=dp)
        end do

        call diagonalize_hermitian_complex_values_only(H, L, eigvals, ierr)

        call check(ierr == 0, "Complex eigenvalues-only should succeed")
        call check(all(eigvals >= -2.0_dp - TOL), "Eigenvalues in band")
        call check(all(eigvals <= 2.0_dp + TOL), "Eigenvalues in band")
    end subroutine test_diag_complex_values_only

end program test_lapack_wrapper
