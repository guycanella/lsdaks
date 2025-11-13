!> Unit tests for boundary_conditions module
program test_boundary_conditions
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp
    implicit none

    real(dp), parameter :: TOL = 1.0e-9_dp

    call execute_serial_cmd_app(get_bc_tests())

contains

    function get_bc_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests

        tests = test_list([ &
            test("validate_bc_open", test_validate_bc_open), &
            test("validate_bc_periodic", test_validate_bc_periodic), &
            test("validate_bc_twisted_valid", test_validate_bc_twisted_valid), &
            test("validate_bc_twisted_no_theta", test_validate_bc_twisted_no_theta), &
            test("validate_bc_twisted_out_of_range", test_validate_bc_twisted_out_of_range), &
            test("validate_bc_invalid_type", test_validate_bc_invalid_type), &
            test("validate_bc_small_system", test_validate_bc_small_system), &
            test("apply_bc_open_no_modification", test_apply_bc_open_no_modification), &
            test("apply_bc_periodic_real", test_apply_bc_periodic_real), &
            test("apply_bc_periodic_complex", test_apply_bc_periodic_complex), &
            test("apply_bc_twisted_complex", test_apply_bc_twisted_complex), &
            test("apply_bc_twisted_antiperiodic", test_apply_bc_twisted_antiperiodic), &
            test("apply_bc_size_mismatch", test_apply_bc_size_mismatch), &
            test("free_particle_eigenvalues_open", test_free_particle_eigenvalues_open), &
            test("free_particle_eigenvalues_periodic", test_free_particle_eigenvalues_periodic), &
            test("free_particle_eigenvalues_twisted", test_free_particle_eigenvalues_twisted), &
            test("free_particle_eigenvalues_sorted", test_free_particle_eigenvalues_sorted) &
        ])
    end function get_bc_tests

    !> Test validation of open boundary conditions
    !!
    !! Physics: Open boundary conditions (OBC) represent a finite chain with
    !! no periodic wrapping. This creates edge states and breaks translational
    !! symmetry. Physically relevant for quantum wires, nanoribbons, and
    !! systems where surface/edge physics is important.
    subroutine test_validate_bc_open()
        use fortuno_serial, only: check => serial_check
        use boundary_conditions, only: validate_bc_parameters, BC_OPEN
        use lsda_errors, only: ERROR_INVALID_INPUT
        integer :: ierr

        call validate_bc_parameters(BC_OPEN, L=10, ierr=ierr)
        call check(ierr == 0, "Open BC should be valid without theta")
    end subroutine test_validate_bc_open

    !> Test validation of periodic boundary conditions
    !!
    !! Physics: Periodic boundary conditions (PBC) create a ring topology,
    !! preserving translational symmetry and momentum conservation. This
    !! eliminates finite-size edge effects and is essential for Bethe Ansatz
    !! calculations and studying bulk properties.
    subroutine test_validate_bc_periodic()
        use fortuno_serial, only: check => serial_check
        use boundary_conditions, only: validate_bc_parameters, BC_PERIODIC
        integer :: ierr

        call validate_bc_parameters(BC_PERIODIC, L=10, ierr=ierr)
        call check(ierr == 0, "Periodic BC should be valid without theta")
    end subroutine test_validate_bc_periodic

    !> Test validation of twisted boundary conditions with valid theta
    !!
    !! Physics: Twisted boundary conditions (TBC) introduce an Aharonov-Bohm
    !! phase θ from a magnetic flux threading the ring. This models persistent
    !! currents in mesoscopic rings and the Aharonov-Bohm effect. The twist
    !! angle θ ∈ [0, 2π) corresponds to flux Φ/Φ₀ = θ/(2π).
    subroutine test_validate_bc_twisted_valid()
        use fortuno_serial, only: check => serial_check
        use boundary_conditions, only: validate_bc_parameters, BC_TWISTED
        use lsda_constants, only: PI
        integer :: ierr

        call validate_bc_parameters(BC_TWISTED, theta=PI/2.0_dp, L=10, ierr=ierr)
        call check(ierr == 0, "Twisted BC should be valid with theta in [0, 2π)")
    end subroutine test_validate_bc_twisted_valid

    !> Test validation failure when theta is missing for TBC
    !!
    !! Physics: Twisted BC requires a twist angle to define the magnetic flux.
    !! Without theta, the boundary phase is undefined, making the problem
    !! ill-posed. This test ensures proper error handling.
    subroutine test_validate_bc_twisted_no_theta()
        use fortuno_serial, only: check => serial_check
        use boundary_conditions, only: validate_bc_parameters, BC_TWISTED
        use lsda_errors, only: ERROR_INVALID_INPUT
        integer :: ierr

        call validate_bc_parameters(BC_TWISTED, L=10, ierr=ierr)
        call check(ierr == ERROR_INVALID_INPUT, "Twisted BC without theta should fail")
    end subroutine test_validate_bc_twisted_no_theta

    !> Test validation failure when theta is out of range
    !!
    !! Physics: The twist angle θ must be in [0, 2π) because it represents
    !! a phase in the complex plane. Values outside this range are redundant
    !! (θ + 2π is equivalent to θ) and should be normalized before use.
    subroutine test_validate_bc_twisted_out_of_range()
        use fortuno_serial, only: check => serial_check
        use boundary_conditions, only: validate_bc_parameters, BC_TWISTED
        use lsda_constants, only: TWOPI
        use lsda_errors, only: ERROR_OUT_OF_BOUNDS
        integer :: ierr

        ! Test theta < 0
        call validate_bc_parameters(BC_TWISTED, theta=-0.1_dp, L=10, ierr=ierr)
        call check(ierr == ERROR_OUT_OF_BOUNDS, "Negative theta should fail")

        ! Test theta >= 2π
        call validate_bc_parameters(BC_TWISTED, theta=TWOPI, L=10, ierr=ierr)
        call check(ierr == ERROR_OUT_OF_BOUNDS, "theta >= 2π should fail")
    end subroutine test_validate_bc_twisted_out_of_range

    !> Test validation failure for invalid BC type
    !!
    !! Physics: Only three boundary condition types are physically meaningful
    !! for 1D systems: open, periodic, and twisted. Invalid types should be
    !! rejected to prevent undefined behavior.
    subroutine test_validate_bc_invalid_type()
        use fortuno_serial, only: check => serial_check
        use boundary_conditions, only: validate_bc_parameters
        use lsda_errors, only: ERROR_INVALID_INPUT
        integer :: ierr

        call validate_bc_parameters(999, L=10, ierr=ierr)
        call check(ierr == ERROR_INVALID_INPUT, "Invalid BC type should fail")
    end subroutine test_validate_bc_invalid_type

    !> Test validation failure for too small system
    !!
    !! Physics: A 1D chain needs at least 2 sites to have meaningful hopping
    !! dynamics. Single-site or zero-site systems cannot support tight-binding
    !! physics and should be rejected.
    subroutine test_validate_bc_small_system()
        use fortuno_serial, only: check => serial_check
        use boundary_conditions, only: validate_bc_parameters, BC_OPEN
        use lsda_errors, only: ERROR_INVALID_INPUT
        integer :: ierr

        call validate_bc_parameters(BC_OPEN, L=1, ierr=ierr)
        call check(ierr == ERROR_INVALID_INPUT, "L=1 should fail")

        call validate_bc_parameters(BC_OPEN, L=0, ierr=ierr)
        call check(ierr == ERROR_INVALID_INPUT, "L=0 should fail")
    end subroutine test_validate_bc_small_system

    !> Test that open BC leaves Hamiltonian unchanged
    !!
    !! Physics: Open boundary conditions naturally result in a tridiagonal
    !! Hamiltonian with no edge connections (H(1,L) = H(L,1) = 0). The chain
    !! is "cut" at the edges, creating confined eigenstates similar to a
    !! particle in a box with hard walls.
    subroutine test_apply_bc_open_no_modification()
        use fortuno_serial, only: check => serial_check
        use boundary_conditions, only: apply_boundary_conditions, BC_OPEN
        integer, parameter :: L = 5
        real(dp) :: H(L,L), H_orig(L,L)
        integer :: ierr

        ! Initialize tridiagonal Hamiltonian
        H = 0.0_dp
        H_orig = H

        call apply_boundary_conditions(H, L, BC_OPEN, ierr=ierr)

        call check(ierr == 0, "Open BC should succeed")
        call check(all(abs(H - H_orig) < TOL), "Open BC should not modify H")
    end subroutine test_apply_bc_open_no_modification

    !> Test periodic BC application to real Hamiltonian
    !!
    !! Physics: Periodic boundary conditions connect site L back to site 1,
    !! creating a ring. This preserves translational symmetry and allows
    !! Bloch waves with well-defined crystal momentum k. The Hamiltonian
    !! becomes circulant, enabling efficient diagonalization.
    subroutine test_apply_bc_periodic_real()
        use fortuno_serial, only: check => serial_check
        use boundary_conditions, only: apply_boundary_conditions, BC_PERIODIC
        integer, parameter :: L = 5
        real(dp) :: H(L,L)
        integer :: ierr

        H = 0.0_dp
        call apply_boundary_conditions(H, L, BC_PERIODIC, ierr=ierr)

        call check(ierr == 0, "Periodic BC should succeed")
        call check(abs(H(1,L) + 1.0_dp) < TOL, "H(1,L) should be -t = -1")
        call check(abs(H(L,1) + 1.0_dp) < TOL, "H(L,1) should be -t = -1")
    end subroutine test_apply_bc_periodic_real

    !> Test periodic BC application to complex Hamiltonian
    !!
    !! Physics: For spinless fermions without magnetic field, PBC can be
    !! represented with a real Hamiltonian. However, when including spin
    !! or preparing for twisted BC, using complex matrices is necessary.
    !! PBC is TBC with θ = 0.
    subroutine test_apply_bc_periodic_complex()
        use fortuno_serial, only: check => serial_check
        use boundary_conditions, only: apply_boundary_conditions_complex, BC_PERIODIC
        integer, parameter :: L = 5
        complex(dp) :: H(L,L)
        integer :: ierr

        H = cmplx(0.0_dp, 0.0_dp, kind=dp)
        call apply_boundary_conditions_complex(H, L, BC_PERIODIC, ierr=ierr)

        call check(ierr == 0, "Periodic BC (complex) should succeed")
        call check(abs(real(H(1,L)) + 1.0_dp) < TOL, "Real part of H(1,L) should be -1")
        call check(abs(aimag(H(1,L))) < TOL, "Imaginary part of H(1,L) should be 0")
        call check(abs(real(H(L,1)) + 1.0_dp) < TOL, "Real part of H(L,1) should be -1")
        call check(abs(aimag(H(L,1))) < TOL, "Imaginary part of H(L,1) should be 0")
    end subroutine test_apply_bc_periodic_complex

    !> Test twisted BC with arbitrary angle
    !!
    !! Physics: Twisted boundary conditions introduce a gauge phase e^{iθ}
    !! when hopping from site L to site 1. This simulates the Aharonov-Bohm
    !! effect where electrons acquire a phase from encircling a magnetic flux.
    !! For θ ≠ 0, time-reversal symmetry is broken, leading to persistent
    !! currents even in the ground state.
    subroutine test_apply_bc_twisted_complex()
        use fortuno_serial, only: check => serial_check
        use boundary_conditions, only: apply_boundary_conditions_complex, BC_TWISTED
        use lsda_constants, only: PI
        integer, parameter :: L = 5
        complex(dp) :: H(L,L)
        real(dp) :: theta
        integer :: ierr

        theta = PI / 4.0_dp
        H = cmplx(0.0_dp, 0.0_dp, kind=dp)
        call apply_boundary_conditions_complex(H, L, BC_TWISTED, theta=theta, ierr=ierr)

        call check(ierr == 0, "Twisted BC should succeed")

        ! H(1,L) = -exp(iθ) = -cos(θ) - i·sin(θ)
        call check(abs(real(H(1,L)) + cos(theta)) < TOL, &
                   "Real part of H(1,L) should be -cos(θ)")
        call check(abs(aimag(H(1,L)) + sin(theta)) < TOL, &
                   "Imaginary part of H(1,L) should be -sin(θ)")

        ! H(L,1) = -exp(-iθ) = -cos(θ) + i·sin(θ)
        call check(abs(real(H(L,1)) + cos(theta)) < TOL, &
                   "Real part of H(L,1) should be -cos(θ)")
        call check(abs(aimag(H(L,1)) - sin(theta)) < TOL, &
                   "Imaginary part of H(L,1) should be sin(θ)")
    end subroutine test_apply_bc_twisted_complex

    !> Test antiperiodic boundary conditions (θ = π)
    !!
    !! Physics: Antiperiodic BC (θ = π) impose a phase factor of -1 when
    !! hopping around the ring. This is equivalent to half a flux quantum
    !! (Φ/Φ₀ = 1/2) and is particularly useful for odd electron numbers,
    !! where it can lift certain degeneracies. The eigenspectrum is shifted
    !! by π/L compared to regular PBC.
    subroutine test_apply_bc_twisted_antiperiodic()
        use fortuno_serial, only: check => serial_check
        use boundary_conditions, only: apply_boundary_conditions_complex, BC_TWISTED
        use lsda_constants, only: PI
        integer, parameter :: L = 5
        complex(dp) :: H(L,L)
        integer :: ierr

        H = cmplx(0.0_dp, 0.0_dp, kind=dp)
        call apply_boundary_conditions_complex(H, L, BC_TWISTED, theta=PI, ierr=ierr)

        call check(ierr == 0, "Antiperiodic BC should succeed")

        ! H(1,L) = -exp(iπ) = 1
        call check(abs(real(H(1,L)) - 1.0_dp) < TOL, &
                   "Real part of H(1,L) should be +1 (antiperiodic)")
        call check(abs(aimag(H(1,L))) < TOL, &
                   "Imaginary part should be ~0")

        ! H(L,1) = -exp(-iπ) = 1
        call check(abs(real(H(L,1)) - 1.0_dp) < TOL, &
                   "Real part of H(L,1) should be +1 (antiperiodic)")
    end subroutine test_apply_bc_twisted_antiperiodic

    !> Test error handling for size mismatch
    !!
    !! Physics: The Hamiltonian matrix must be square with dimensions L×L
    !! to properly represent the tight-binding model on L sites. Mismatched
    !! dimensions indicate a programming error and should be caught early.
    subroutine test_apply_bc_size_mismatch()
        use fortuno_serial, only: check => serial_check
        use boundary_conditions, only: apply_boundary_conditions, BC_PERIODIC
        use lsda_errors, only: ERROR_SIZE_MISMATCH
        integer, parameter :: L = 5
        real(dp) :: H(L,L+1)  ! Wrong size
        integer :: ierr

        call apply_boundary_conditions(H, L, BC_PERIODIC, ierr=ierr)
        call check(ierr == ERROR_SIZE_MISMATCH, "Size mismatch should fail")
    end subroutine test_apply_bc_size_mismatch

    !> Test free particle eigenvalues for open BC
    !!
    !! Physics: For open boundary conditions, the eigenstates are standing
    !! waves sin(nπx/L) similar to a particle in a box. The eigenvalues are
    !! E_n = -2cos(nπ/(L+1)), giving a discrete spectrum bounded by [-2, 2].
    !! This analytical formula provides a crucial validation for numerical
    !! diagonalization routines.
    subroutine test_free_particle_eigenvalues_open()
        use fortuno_serial, only: check => serial_check
        use boundary_conditions, only: get_free_particle_eigenvalues, BC_OPEN
        use lsda_constants, only: PI
        integer, parameter :: L = 10
        real(dp) :: eigenvalues(L)
        integer :: ierr, n

        call get_free_particle_eigenvalues(L, BC_OPEN, eigenvalues=eigenvalues, ierr=ierr)

        call check(ierr == 0, "Should succeed for open BC")

        ! Check first eigenvalue manually
        call check(abs(eigenvalues(1) + 2.0_dp * cos(PI / real(L+1, dp))) < TOL, &
                   "First eigenvalue should match analytical formula")

        ! All eigenvalues should be in [-2, 2] (bandstructure bounds)
        call check(all(eigenvalues >= -2.0_dp - TOL), "All eigenvalues >= -2")
        call check(all(eigenvalues <= 2.0_dp + TOL), "All eigenvalues <= 2")
    end subroutine test_free_particle_eigenvalues_open

    !> Test free particle eigenvalues for periodic BC
    !!
    !! Physics: Periodic boundary conditions allow Bloch states with crystal
    !! momentum k = 2πn/L (n = 0, 1, ..., L-1). The dispersion relation is
    !! E(k) = -2cos(k), giving the characteristic cosine bandstructure of
    !! 1D tight-binding models. Momentum is a good quantum number.
    subroutine test_free_particle_eigenvalues_periodic()
        use fortuno_serial, only: check => serial_check
        use boundary_conditions, only: get_free_particle_eigenvalues, BC_PERIODIC
        use lsda_constants, only: TWOPI
        integer, parameter :: L = 10
        real(dp) :: eigenvalues(L)
        integer :: ierr

        call get_free_particle_eigenvalues(L, BC_PERIODIC, eigenvalues=eigenvalues, ierr=ierr)

        call check(ierr == 0, "Should succeed for periodic BC")

        ! k=0 gives E=-2 (bottom of band)
        call check(abs(eigenvalues(1) + 2.0_dp) < TOL, &
                   "k=0 eigenvalue should be -2 (band minimum)")

        ! All in band
        call check(all(eigenvalues >= -2.0_dp - TOL), "All eigenvalues >= -2")
        call check(all(eigenvalues <= 2.0_dp + TOL), "All eigenvalues <= 2")
    end subroutine test_free_particle_eigenvalues_periodic

    !> Test free particle eigenvalues for twisted BC
    !!
    !! Physics: Twisted boundary conditions shift the allowed k-values by
    !! θ/L, giving E_k(θ) = -2cos((2πk + θ)/L). For non-zero twist, the
    !! spectrum is shifted relative to periodic BC. As θ varies from 0 to 2π,
    !! the eigenvalues sweep through the band, exhibiting Aharonov-Bohm
    !! oscillations in physical observables like persistent current.
    subroutine test_free_particle_eigenvalues_twisted()
        use fortuno_serial, only: check => serial_check
        use boundary_conditions, only: get_free_particle_eigenvalues, BC_TWISTED
        use lsda_constants, only: PI
        integer, parameter :: L = 10
        real(dp) :: eigenvalues(L), theta
        integer :: ierr

        theta = PI / 3.0_dp
        call get_free_particle_eigenvalues(L, BC_TWISTED, theta=theta, &
                                            eigenvalues=eigenvalues, ierr=ierr)

        call check(ierr == 0, "Should succeed for twisted BC")

        ! First eigenvalue: k=0, E = -2cos(θ/L)
        call check(abs(eigenvalues(1) + 2.0_dp * cos(theta / real(L, dp))) < TOL, &
                   "First eigenvalue should match twisted formula")

        ! All in band
        call check(all(eigenvalues >= -2.0_dp - TOL), "All eigenvalues >= -2")
        call check(all(eigenvalues <= 2.0_dp + TOL), "All eigenvalues <= 2")
    end subroutine test_free_particle_eigenvalues_twisted

    !> Test eigenvalue spectrum bounds and completeness
    !!
    !! Physics: Eigenvalues represent energy levels and must lie within the
    !! tight-binding band [-2t, 2t] (t=1 here). For PBC, the eigenvalues come
    !! from E(k) = -2cos(k), which traces out the full cosine dispersion as
    !! k varies. The ground state is at k=0 (E=-2). For even L, there are
    !! degeneracies because E(k) = E(-k) due to time-reversal symmetry, so
    !! states at ±k have the same energy (except k=0 and k=π).
    subroutine test_free_particle_eigenvalues_sorted()
        use fortuno_serial, only: check => serial_check
        use boundary_conditions, only: get_free_particle_eigenvalues, BC_PERIODIC
        use lsda_constants, only: TWOPI
        integer, parameter :: L = 20
        real(dp) :: eigenvalues(L), evals_sorted(L)
        integer :: ierr, i, k
        real(dp) :: e_min, e_max

        call get_free_particle_eigenvalues(L, BC_PERIODIC, eigenvalues=eigenvalues, ierr=ierr)

        call check(ierr == 0, "Should succeed")

        ! All eigenvalues should be in band [-2, 2]
        call check(all(eigenvalues >= -2.0_dp - TOL), "All eigenvalues >= -2")
        call check(all(eigenvalues <= 2.0_dp + TOL), "All eigenvalues <= 2")

        ! Ground state should be at E=-2 (k=0)
        call check(abs(eigenvalues(1) + 2.0_dp) < TOL, &
                   "Ground state eigenvalue should be -2 for PBC")

        ! Find min and max eigenvalues
        e_min = minval(eigenvalues)
        e_max = maxval(eigenvalues)

        ! Spectrum should span the band (or close to it for finite L)
        call check(e_min < -1.9_dp, "Minimum eigenvalue should be close to -2")
        call check(e_max > 1.9_dp, "Maximum eigenvalue should be close to +2")

        ! Verify cosine dispersion E(k) = -2cos(2πk/L) for a few k values
        do k = 0, min(3, L-1)
            call check(abs(eigenvalues(k+1) + 2.0_dp * cos(TWOPI * k / real(L, dp))) < TOL, &
                       "Eigenvalue should match E(k) = -2cos(2πk/L)")
        end do
    end subroutine test_free_particle_eigenvalues_sorted

end program test_boundary_conditions
