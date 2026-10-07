!> Release benchmark for the partial open-chain diagonalization path.
!!
!! Compare the production DSTEVR path against the legacy full-spectrum dense
!! wrapper at L = 100, 500 and 1000. Timing is machine dependent, so this is
!! an executable rather than a CI-gating unit test.
program benchmark_partial_diagonalization
    use lsda_constants, only: dp
    use lsda_errors, only: ERROR_SUCCESS
    use lapack_wrapper, only: diag_workspace_t, cleanup_diag_workspace, &
                              diagonalize_open_tridiagonal, diagonalize_symmetric_real
    implicit none

    integer, parameter :: L_LIST(3) = [100, 500, 1000]
    integer :: i, j, L, n_vec, ierr
    real(dp) :: t0, t1, partial_seconds, full_seconds, max_error
    real(dp), allocatable :: potential(:), hamiltonian(:,:), partial_values(:), full_values(:)
    real(dp), allocatable :: partial_vectors(:,:), full_vectors(:,:)
    type(diag_workspace_t) :: workspace

    print '(A)', ' L       partial(s)      full(s)       speedup     max |delta E|'
    do i = 1, size(L_LIST)
        L = L_LIST(i)
        ! Low filling is the target workload: the SCF needs the occupied
        ! states plus a small Fermi-shell buffer, not a fixed fraction of L.
        n_vec = min(L, 55)
        allocate(potential(L), hamiltonian(L, L), partial_values(n_vec), full_values(L), &
                 partial_vectors(L, n_vec), full_vectors(L, L))

        do j = 1, L
            potential(j) = 0.1_dp * sin(real(j, dp))
        end do
        hamiltonian = 0.0_dp
        do j = 1, L
            hamiltonian(j, j) = potential(j)
        end do
        do j = 1, L - 1
            hamiltonian(j, j + 1) = -1.0_dp
            hamiltonian(j + 1, j) = -1.0_dp
        end do

        call cpu_time(t0)
        call diagonalize_open_tridiagonal(potential, L, n_vec, partial_values, partial_vectors, workspace, ierr)
        call cpu_time(t1)
        partial_seconds = t1 - t0
        if (ierr /= ERROR_SUCCESS) error stop 'partial diagonalization failed'

        call cpu_time(t0)
        call diagonalize_symmetric_real(hamiltonian, L, full_values, full_vectors, ierr)
        call cpu_time(t1)
        full_seconds = t1 - t0
        if (ierr /= ERROR_SUCCESS) error stop 'full diagonalization failed'

        max_error = maxval(abs(partial_values - full_values(1:n_vec)))
        print '(I5,3X,F10.4,3X,F10.4,3X,F9.2,3X,ES12.4)', L, partial_seconds, full_seconds, &
            full_seconds / max(partial_seconds, tiny(1.0_dp)), max_error

        call cleanup_diag_workspace(workspace)
        deallocate(potential, hamiltonian, partial_values, full_values, partial_vectors, full_vectors)
    end do
end program benchmark_partial_diagonalization
